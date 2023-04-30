#include "renderer.h"

#include <algorithm> //sort

#include "camera.h"
#include "../gfx/gfx.h"
#include "../gfx/shader.h"
#include "../gfx/mesh.h"
#include "../gfx/texture.h"
#include "../gfx/fbo.h"
#include "../pipeline/prefab.h"
#include "../pipeline/material.h"
#include "../pipeline/animation.h"
#include "../utils/utils.h"
#include "../extra/hdre.h"
#include "../core/ui.h"

#include "scene.h"


using namespace SCN;

//some globals
GFX::Mesh sphere;

Renderer::Renderer(const char* shader_atlas_filename)
{
	render_wireframe = false;
	render_boundaries = false;
	scene = nullptr;
	skybox_cubemap = nullptr;
	render_mode = eRenderMode::LIGHTSINGLEPASS;

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
}

void Renderer::setupScene()
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;

	lights.clear();

	//render entities
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible)
			continue;
		if (ent->getType() == eEntityType::PREFAB) {
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
		}
		else if (ent->getType() == eEntityType::LIGHT) {
			lights.push_back((SCN::LightEntity*)ent);
		}
	}

}

float SCN::Renderer::distance(vec3 p1, vec3 p2)
{
	return std::sqrt(pow(p1.x-p2.x,2.0)+ pow(p1.y - p2.y, 2.0)+ pow(p1.z - p2.z, 2.0));
}



void Renderer::renderScene(SCN::Scene* scene, Camera* camera)
{
	this->scene = scene;
	setupScene();

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if(skybox_cubemap)
		renderSkybox(skybox_cubemap);

	//Restart the list that stores the elements in our scene so that we don't store all elements in every render iteration
	renderCalls.clear();

	//render entities
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible )
			continue;

		//is a prefab!
		if (ent->getType() == eEntityType::PREFAB)
		{
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
			if (pent->prefab) {

				renderNode( &pent->root, camera);
				
			}

		}
	}

	std::sort(renderCalls.begin(), renderCalls.end(), [](const RenderCall& a, const RenderCall& b) {
		return a.distance_to_camera > b.distance_to_camera;
	});

	switch (render_mode)
	{
	case SCN::FLAT:
		for (const RenderCall& obj : renderCalls) {
			renderMeshWithMaterial(obj.model, obj.mesh, obj.material);
		}
		break;
	case SCN::LIGHTMULTIPASS:
		for (const RenderCall& obj : renderCalls) {
			renderMeshWithMaterialLight(obj.model, obj.mesh, obj.material);
		}
	case SCN::LIGHTSINGLEPASS:
		for (const RenderCall& obj : renderCalls) {
			renderMeshWithMaterialLightSinglePass(obj.model, obj.mesh, obj.material);
		}
	default:
		break;
	}
	
	

	


}


void Renderer::renderSkybox(GFX::Texture* cubemap)
{
	Camera* camera = Camera::current;

	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_CULL_FACE);
	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	GFX::Shader* shader = GFX::Shader::Get("skybox");
	if (!shader)
		return;
	shader->enable();

	Matrix44 m;
	m.setTranslation(camera->eye.x, camera->eye.y, camera->eye.z);
	m.scale(10, 10, 10);
	shader->setUniform("u_model", m);
	cameraToShader(camera, shader);
	shader->setUniform("u_texture", cubemap, 0);
	sphere.render(GL_TRIANGLES);
	shader->disable();
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glEnable(GL_DEPTH_TEST);
}

//renders a node of the prefab and its children
void Renderer::renderNode(SCN::Node* node, Camera* camera)
{
	if (!node->visible)
		return;

	//compute global matrix
	Matrix44 node_model = node->getGlobalMatrix(true);

	//does this node have a mesh? then we must render it
	if (node->mesh && node->material)
	{
		//compute the bounding box of the object in world space (by using the mesh bounding box transformed to world space)
		BoundingBox world_bounding = transformBoundingBox(node_model,node->mesh->box);
		
		//if bounding box is inside the camera frustum then the object is probably visible
		if (camera->testBoxInFrustum(world_bounding.center, world_bounding.halfsize) )
		{
			if(render_boundaries)
				node->mesh->renderBounding(node_model, true);

			RenderCall rc;
			
			vec3 nodepos = node->getGlobalMatrix().getTranslation();
			rc.mesh = node->mesh;
			rc.material = node->material;
			rc.model = node_model;
			rc.distance_to_camera = distance(camera->eye, nodepos);
			

			renderCalls.push_back(rc);

			// We don't want to render the object yet, we will store each object and its distance to camera
			// so that later we can sort from furthest to closest distances.
			//renderMeshWithMaterial(node_model, node->mesh, node->material);
		}
	}

	//iterate recursively with children
	for (int i = 0; i < node->children.size(); ++i)
		renderNode( node->children[i], camera);
}

//renders a mesh given its transform and material
void Renderer::renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material )
		return;
    assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	//GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;
	

	GFX::Texture*  white_texture = GFX::Texture::getWhiteTexture();
	GFX::Texture* albedo_texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	//texture = material->metallic_roughness_texture;
	//texture = material->normal_texture;
	//texture = material->occlusion_texture;
	if (albedo_texture == NULL)
		albedo_texture = white_texture; //a 1x1 white texture

	if (emissive_texture == NULL)
		emissive_texture = white_texture;

	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if(material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
    assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("texture");

    assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t );

	shader->setUniform("u_color", material->color);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_texture", albedo_texture?albedo_texture:white_texture, 0);
	//shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	//shader->setUniform("u_ambient_light",scene->ambient_light);

	if (render_wireframe)
		glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}



//renders a mesh given its transform and material and takes into consideration all the lights in the scene
void Renderer::renderMeshWithMaterialLight(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	//GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;


	GFX::Texture* white_texture = GFX::Texture::getWhiteTexture();
	GFX::Texture* albedo_texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	//texture = material->metallic_roughness_texture;
	//texture = material->normal_texture;
	//texture = material->occlusion_texture;
	if (albedo_texture == NULL)
		albedo_texture = white_texture; //a 1x1 white texture

	if (emissive_texture == NULL)
		emissive_texture = white_texture;

	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("light");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white_texture, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);



	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	shader->setUniform("u_ambient_light", scene->ambient_light);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	glDepthFunc(GL_LEQUAL);

	if (lights.size() == 0) {
		//Take care of scenario
		shader->setUniform("u_light_type", 0);
		mesh->render(GL_TRIANGLES);

	}

	for (int i = 0; i < lights.size(); ++i) {

		LightEntity* light = lights[i];
		shader->setUniform("u_light_position", light->root.model.getTranslation());
		shader->setUniform("u_light_color", light->color * light->intensity);
		shader->setUniform("u_light_info", vec4((int)light->light_type,light->near_distance,light->max_distance,0));
		shader->setUniform("u_light_front", light->root.model.rotateVector(vec3(0,0,1)));
		if(light->light_type == eLightType::SPOT)
			shader->setUniform("u_light_cone", vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD)));

		//if (light->light_type != eLightType::POINT && light->light_type != eLightType::DIRECTIONAL)
		//	continue;
		//do the draw call that renders the mesh into the screen
		mesh->render(GL_TRIANGLES);

		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE);


		shader->setUniform("u_ambient_light", vec3(0.0));
		shader->setUniform("u_emissive_factor", vec3(0.0));

		
	}


	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glDepthFunc(GL_LESS);
}

void Renderer::renderMeshWithMaterialLightSinglePass(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	//GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;


	GFX::Texture* white_texture = GFX::Texture::getWhiteTexture();
	GFX::Texture* albedo_texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	//texture = material->metallic_roughness_texture;
	//texture = material->normal_texture;
	//texture = material->occlusion_texture;
	if (albedo_texture == NULL)
		albedo_texture = white_texture; //a 1x1 white texture

	if (emissive_texture == NULL)
		emissive_texture = white_texture;

	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
	{
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	}
	else
		glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("lightSinglePass");

	assert(glGetError() == GL_NO_ERROR);

	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();

	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);
	float t = getTime();
	shader->setUniform("u_time", t);

	shader->setUniform("u_color", material->color);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white_texture, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);

	//this is used to say which is the alpha threshold to what we should not paint a pixel on the screen (to cut polygons according to texture alpha)
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);

	shader->setUniform("u_ambient_light", scene->ambient_light);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	
	vec3 light_position[4] = {};
	vec3 light_color[4] = {};
	vec4 light_info[4] = {};
	vec3 light_front[4] = {};
	vec2 light_cone[4] = {};

	for (int i = 0; i < lights.size(); ++i) {
		if (i < 4) 
		{
			light_color[i] = lights[i]->color;
			light_position[i] = lights[i]->root.model.getTranslation();
			light_front[i] = lights[i]->root.model.rotateVector(vec3(0, 0, 1));
			light_info[i] = vec4((int)lights[i]->light_type, lights[i]->near_distance, lights[i]->max_distance, 0);//(light_type, near_distance, max_distance, xx)
			if (lights[i]->light_type == eLightType::SPOT)
				light_cone[i] = vec2(cos(lights[i]->cone_info.x * DEG2RAD), cos(lights[i]->cone_info.y * DEG2RAD));

		}
			
	}

	shader->setUniform("u_num_lights", (int)lights.size());
	shader->setUniform3Array("u_light_position", (float*)&light_position,4);
	shader->setUniform3Array("u_light_color", (float*)&light_color, 4);
	shader->setUniform3Array("u_light_front", (float*)&light_front, 4);
	shader->setUniform4Array("u_light_info", (float*)&light_info, 4);
	shader->setUniform2Array("u_light_cone", (float*)&light_cone, 4);

	
	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);

	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);

	

	//disable shader
	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

}


void SCN::Renderer::cameraToShader(Camera* camera, GFX::Shader* shader)
{
	shader->setUniform("u_viewprojection", camera->viewprojection_matrix );
	shader->setUniform("u_camera_position", camera->eye);
}

#ifndef SKIP_IMGUI

void Renderer::showUI()
{
		
	ImGui::Checkbox("Wireframe", &render_wireframe);
	ImGui::Checkbox("Boundaries", &render_boundaries);
	ImGui::Combo("Render Mode", (int*) & render_mode, "FLAT\0LIGHTSMULTIPASS\0LIGHTSSINGLEPASS\0", 2);

	//add here your stuff
	//...
}

#else
void Renderer::showUI() {}
#endif