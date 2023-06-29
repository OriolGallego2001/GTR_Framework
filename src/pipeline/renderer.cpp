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

GFX::FBO* illumination_fbo = nullptr;

Renderer::Renderer(const char* shader_atlas_filename)
{
	show_volumetric = false;
	air_density = 0.0001;
	render_wireframe = false;
	render_boundaries = false;
	show_gbuffers = false;
	scene = nullptr;
	skybox_cubemap = nullptr;
	render_mode = eRenderMode::DEFERRED;
	show_shadowmap = false;
	show_shadows = true;
	normal_maps = true;
	gbuffers_fbo = nullptr;
	shininess_coef = 10;
	shadow_atlas_fbo = new GFX::FBO();
	shadow_atlas_fbo->setDepthOnly(1024* max_shadowmaps_pow2, 1024* max_shadowmaps_pow2);
	current_shadow_texture = shadow_atlas_fbo->depth_texture;

	tonemapper_scale = 1.0; //color scale before tonemapper
	average_lum = 1.0;
	lumwhite2 = 1.0;
	igamma = 1.0;

	vec2 size = CORE::getWindowSize();
	light_fbo = new GFX::FBO();
	light_fbo->create(size.x, size.y, 1, GL_RGB, GL_HALF_FLOAT);

	ssao_fbo = new GFX::FBO();
	ssao_fbo->create(size.x, size.y, 1, GL_RGB, GL_UNSIGNED_BYTE, false);

	irr_fbo = new GFX::FBO();
	irr_fbo->create(64, 64, 1, GL_RGB, GL_FLOAT, false);

	if (!GFX::Shader::LoadAtlas(shader_atlas_filename))
		exit(1);
	GFX::checkGLErrors();

	sphere.createSphere(1.0f);
	sphere.uploadToVRAM();


	irradiance_multiplier = 1.0;
	show_probes = false;
	probes_texture = nullptr;
}

void Renderer::setupScene(Camera* camera)
{
	if (scene->skybox_filename.size())
		skybox_cubemap = GFX::Texture::Get(std::string(scene->base_folder + "/" + scene->skybox_filename).c_str());
	else
		skybox_cubemap = nullptr;

	//Restart the list that stores the elements in our scene so that we don't store all elements in every render iteration
	opaqueRenderCalls.clear();
	transparentRenderCalls.clear();
	lights.clear();

	//render entities
	for (int i = 0; i < scene->entities.size(); ++i)
	{
		BaseEntity* ent = scene->entities[i];
		if (!ent->visible)
			continue;
		if (ent->getType() == eEntityType::PREFAB) {
			PrefabEntity* pent = (SCN::PrefabEntity*)ent;
			if (pent->prefab) {
				// Camera* camera = Camera::current;
				renderNode(&pent->root, camera);
			}
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
	setupScene(camera);

	if (show_shadows) {
		generateShadowMaps();
		for (auto light : lights) {
			if (light->shadowmap_fbo) {
				light->shadowmap = light->shadowmap_fbo->depth_texture;
			}
		}
	}

	
	std::sort(opaqueRenderCalls.begin(), opaqueRenderCalls.end(), [](const RenderCall& a, const RenderCall& b) {
		return a.distance_to_camera > b.distance_to_camera;
		});
	std::sort(transparentRenderCalls.begin(), transparentRenderCalls.end(), [](const RenderCall& a, const RenderCall& b) {
		return a.distance_to_camera > b.distance_to_camera;
		});

	renderFrame(scene, camera);

	//Debug
	if (show_shadowmap)
		debug_shadowmap();


}
void Renderer::renderFrame(SCN::Scene* scene, Camera* camera) 
{
	
	if (render_mode == eRenderMode::DEFERRED) {
		renderDeferred(scene, camera);
	}
	else {
		renderForward(scene, camera);
	}
	
}
void Renderer::renderForward(SCN::Scene* scene, Camera* camera)
{

	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	camera->enable();

	//set the clear color (the background color)
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);

	// Clear the color and the depth buffer
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	GFX::checkGLErrors();

	//render skybox
	if (skybox_cubemap && render_mode == eRenderMode::FLAT)
		renderSkybox(skybox_cubemap);

	renderRenderCalls(scene, camera);
}
void Renderer::renderRenderCalls(SCN::Scene * scene, Camera * camera)
	{
	

	switch (render_mode)
	{
	case SCN::FLAT:
		for (const RenderCall& obj : opaqueRenderCalls) {
			renderMeshWithMaterialFlat(obj.model, obj.mesh, obj.material);
		}
		for (const RenderCall& obj : transparentRenderCalls) {
			renderMeshWithMaterialFlat(obj.model, obj.mesh, obj.material);
		}
		break;
	case SCN::TEXTURED:
		for (const RenderCall& obj : opaqueRenderCalls) {
			renderMeshWithMaterial(obj.model, obj.mesh, obj.material);
		}
		for (const RenderCall& obj : transparentRenderCalls) {
			renderMeshWithMaterial(obj.model, obj.mesh, obj.material);
		}
		break;
	case SCN::LIGHTMULTIPASS:
		for (const RenderCall& obj : opaqueRenderCalls) {
			renderMeshWithMaterialLight(obj.model, obj.mesh, obj.material);
		}
		for (const RenderCall& obj : transparentRenderCalls) {
			renderMeshWithMaterialLight(obj.model, obj.mesh, obj.material);
		}
		break;
	case SCN::LIGHTSINGLEPASS:
		for (const RenderCall& obj : opaqueRenderCalls) {
			renderMeshWithMaterialLightSinglePass(obj.model, obj.mesh, obj.material);
		}
		for (const RenderCall& obj : transparentRenderCalls) {
			renderMeshWithMaterialLightSinglePass(obj.model, obj.mesh, obj.material);
		}
		break;
	case SCN::DEFERRED:
		for (const RenderCall& obj : opaqueRenderCalls) {
			renderMeshWithMaterialGBuffers(obj.model, obj.mesh, obj.material);
		}
		for (const RenderCall& obj : transparentRenderCalls) {
			renderMeshWithMaterialLight(obj.model, obj.mesh, obj.material);
		}
		break;
	default:
		break;
	}
	
	

	


}
void SCN::Renderer::lightToShader(LightEntity* light, GFX::Shader* shader)
{
	shader->setUniform("u_light_position", light->root.model.getTranslation());
	shader->setUniform("u_light_color", light->color * light->intensity);
	shader->setUniform("u_light_info", vec4((int)light->light_type, light->near_distance, light->max_distance, shininess_coef));
	shader->setUniform("u_light_front", light->root.model.rotateVector(vec3(0, 0, 1)));
	if (light->light_type == eLightType::SPOT)
		shader->setUniform("u_light_cone", vec2(cos(light->cone_info.x * DEG2RAD), cos(light->cone_info.y * DEG2RAD)));


	shader->setUniform("u_shadow_info", vec2(light->shadowmap ? 1 : 0, light->shadow_bias));
	if (light->shadowmap)
	{
		shader->setUniform("u_shadowmap", light->shadowmap, 8);
		shader->setUniform("u_shadow_viewproj", light->shadow_viewproj);

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


			if (node->material->alpha_mode == NO_ALPHA) {
				opaqueRenderCalls.push_back(rc);
			}
			else {
				transparentRenderCalls.push_back(rc);
			}

			

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



//renders a mesh given its transform and material
void Renderer::renderMeshWithMaterialFlat(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);
	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	//GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;
	//select the blending
	if (material->alpha_mode == SCN::eAlphaMode::BLEND)
		return;
	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);
	glEnable(GL_DEPTH_TEST);
	//chose a shader
	shader = GFX::Shader::Get("flat");
	assert(glGetError() == GL_NO_ERROR);
	//no shader? then nothing to render
	if (!shader)
		return;
	shader->enable();
	//upload uniforms
	shader->setUniform("u_model", model);
	cameraToShader(camera, shader);	
	//do the draw call that renders the mesh into the screen
	mesh->render(GL_TRIANGLES);
	//disable shader
	shader->disable();
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
	GFX::Texture* metalic_texture = material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	GFX::Texture* normalmap_texture = material->textures[SCN::eTextureChannel::NORMALMAP].texture;
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
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white_texture, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);
	shader->setUniform("u_occlusion_texture", metalic_texture ? metalic_texture : white_texture, 2);
	shader->setUniform("u_normalmap_texture", normal_maps?normalmap_texture : white_texture);
	shader->setUniform("u_roughness", material->roughness_factor);
	shader->setUniform("u_metalness", material->metallic_factor);


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
		lightToShader(light, shader);
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
	GFX::Texture* normalmap_texture = material->textures[SCN::eTextureChannel::NORMALMAP].texture;
	GFX::Texture* metalic_texture = material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	

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

	glDepthFunc(GL_LEQUAL);

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
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white_texture, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);
	shader->setUniform("u_occlusion_texture", metalic_texture ? metalic_texture : white_texture, 2);
	shader->setUniform("u_normalmap_texture", normal_maps ? normalmap_texture : white_texture);

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
	vec3 light_atlas_info[4] = {};
	mat4 light_shadow_viewproj[4] = {};

	for (int i = 0; i < lights.size(); ++i) {
		if (i < 4) 
		{
			light_color[i] = lights[i]->color;
			light_position[i] = lights[i]->root.model.getTranslation();
			light_front[i] = lights[i]->root.model.rotateVector(vec3(0, 0, 1));
			light_info[i] = vec4((int)lights[i]->light_type, lights[i]->near_distance, lights[i]->max_distance, shininess_coef);//(light_type, near_distance, max_distance, xx)
			if (lights[i]->light_type == eLightType::SPOT)
				light_cone[i] = vec2(cos(lights[i]->cone_info.x * DEG2RAD), cos(lights[i]->cone_info.y * DEG2RAD));

			light_atlas_info[i] = vec3(1024 * i, 1024 * i,lights[i]->shadow_bias);
			light_shadow_viewproj[i] = lights[i]->shadow_viewproj;

		}
			
	}

	shader->setUniform("u_num_lights", (int)lights.size());
	shader->setUniform3Array("u_light_position", (float*)&light_position,4);
	shader->setUniform3Array("u_light_color", (float*)&light_color, 4);
	shader->setUniform3Array("u_light_front", (float*)&light_front, 4);
	shader->setUniform4Array("u_light_info", (float*)&light_info, 4);
	shader->setUniform2Array("u_light_cone", (float*)&light_cone, 4);
	
	shader->setUniform("u_shadow_atlas", current_shadow_texture);
	shader->setUniform3Array("u_shadow_info", (float*)&light_atlas_info, 4);
	shader->setMatrix44Array("u_shadow_viewproj", light_shadow_viewproj, 4);
	
	

	
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

void SCN::Renderer::renderMeshWithMaterialGBuffers(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material)
{
	//in case there is nothing to do
	if (!mesh || !mesh->getNumVertices() || !material)
		return;
	assert(glGetError() == GL_NO_ERROR);

	if (material->alpha_mode == eAlphaMode::BLEND)
		return;

	//define locals to simplify coding
	GFX::Shader* shader = NULL;
	//GFX::Texture* texture = NULL;
	Camera* camera = Camera::current;


	GFX::Texture* white_texture = GFX::Texture::getWhiteTexture();
	GFX::Texture* albedo_texture = material->textures[SCN::eTextureChannel::ALBEDO].texture;
	GFX::Texture* emissive_texture = material->textures[SCN::eTextureChannel::EMISSIVE].texture;
	GFX::Texture* metalic_texture = material->textures[SCN::eTextureChannel::METALLIC_ROUGHNESS].texture;
	GFX::Texture* normalmap_texture = material->textures[SCN::eTextureChannel::NORMALMAP].texture;

	//texture = material->metallic_roughness_texture;
	//texture = material->normal_texture;
	//texture = material->occlusion_texture;
	if (albedo_texture == NULL)
		albedo_texture = white_texture; //a 1x1 white texture

	if (emissive_texture == NULL)
		emissive_texture = white_texture;

	
	glDisable(GL_BLEND);

	//select if render both sides of the triangles
	if (material->two_sided)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);
	assert(glGetError() == GL_NO_ERROR);

	glEnable(GL_DEPTH_TEST);

	//chose a shader
	shader = GFX::Shader::Get("gbuffers");

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
	//shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_emissive_factor", material->emissive_factor);
	shader->setUniform("u_albedo_texture", albedo_texture ? albedo_texture : white_texture, 0);
	shader->setUniform("u_emissive_texture", emissive_texture ? emissive_texture : white_texture, 1);
	shader->setUniform("u_roughness", material->roughness_factor);
	shader->setUniform("u_metalness", material->metallic_factor);

	//shader->setUniform("u_occlusion_texture", metalic_texture ? metalic_texture : white_texture, 2);
	//shader->setUniform("u_normalmap_texture", normal_maps ? normalmap_texture : white_texture);
	shader->setUniform("u_alpha_cutoff", material->alpha_mode == SCN::eAlphaMode::MASK ? material->alpha_cutoff : 0.001f);
	shader->setUniform("u_metalness_texture", metalic_texture);

	//shader->setUniform("u_ambient_light", scene->ambient_light);

	if (render_wireframe)
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);


	mesh->render(GL_TRIANGLES);

	shader->disable();

	//set the render state as it was before to avoid problems with future renders
	glDisable(GL_BLEND);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}


void Renderer::renderDeferred(Scene* scene, Camera* camera)
{
	vec2 size = CORE::getWindowSize();
	GFX::Mesh* quad = GFX::Mesh::getQuad();
	GFX::Shader* shader = nullptr;

	if (!gbuffers_fbo)
	{
		gbuffers_fbo = new GFX::FBO();
		gbuffers_fbo->create(size.x, size.y, 3, GL_RGBA, GL_UNSIGNED_BYTE, true);

		volumetric_fbo = new GFX::FBO();
		volumetric_fbo->create(size.x, size.y, 1, GL_RGBA, GL_UNSIGNED_BYTE,false);
	}
	gbuffers_fbo->bind();


	//gbuffers_fbo->enableBuffers(true, false, false, false);
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	//gbuffers_fbo->enableAllBuffers();
	camera->enable();
	renderRenderCalls(scene, camera);
	gbuffers_fbo->unbind();



	volumetric_fbo->bind();


	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);


	shader = GFX::Shader::Get("volumetric");
	shader->enable();
	

	shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);
	shader->setUniform("u_camera_position", camera->eye);
	shader->setUniform("u_iRes", vec2(1 / size.x, 1 / size.y));
	shader->setUniform("u_ivp", camera->inverse_viewprojection_matrix);
	shader->setUniform("u_air_density", air_density);

	LightEntity* spot = lights[0];
	
	lightToShader(spot, shader);
	quad->render(GL_TRIANGLES);

	

	volumetric_fbo->unbind();

	light_fbo->bind();
	camera->enable();
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);
	glClearColor(scene->background_color.x, scene->background_color.y, scene->background_color.z, 1.0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	if(skybox_cubemap)
		renderSkybox(skybox_cubemap);

	glDisable(GL_DEPTH_TEST);
	shader = GFX::Shader::Get("deferred_global");
	shader->enable();
	
	shader->setTexture("u_albedo_texture",gbuffers_fbo->color_textures[0],0);
	shader->setTexture("u_normal_texture", gbuffers_fbo->color_textures[1], 1);
	shader->setTexture("u_extra_texture", gbuffers_fbo->color_textures[2], 2);
	shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);

	shader->setUniform("u_ambient_light", scene->ambient_light);

	quad->render(GL_TRIANGLES);

	shader->disable();

	shader = GFX::Shader::Get("deferred_light");
	shader->enable();

	shader->setTexture("u_albedo_texture", gbuffers_fbo->color_textures[0], 0);
	shader->setTexture("u_normal_texture", gbuffers_fbo->color_textures[1], 1);
	shader->setTexture("u_extra_texture", gbuffers_fbo->color_textures[2], 2);
	shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);

	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	shader->setUniform("u_iRes", vec2(1 / size.x, 1 / size.y));
	shader->setUniform("u_inv_viewproj", camera->inverse_viewprojection_matrix);
	for (auto light : lights)
	{	
		if (light->light_type == eLightType::POINT || light->light_type == eLightType::SPOT)
		{
			mat4 m;
			vec3 lightpos = light->root.model.getTranslation();
			m.setTranslation(lightpos.x, lightpos.y, lightpos.z);
			m.scale(light->max_distance, light->max_distance, light->max_distance);
			shader->setUniform("u_model", m);
			cameraToShader(camera, shader);
			glFrontFace(GL_CW);
			lightToShader(light, shader);
			
			sphere.render(GL_TRIANGLES);
			glFrontFace(GL_CCW);
		}
		else {
			lightToShader(light, shader);
			shader->setUniform("u_camera_position", camera->eye);
			shader->setUniform("u_viewprojection", mat4::IDENTITY);
			shader->setUniform("u_model", mat4::IDENTITY);
			quad->render(GL_TRIANGLES);
		}
	}
	glDisable(GL_BLEND);

	if (volumetric_fbo)
	{
		glDisable(GL_DEPTH_TEST);
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		volumetric_fbo->color_textures[0]->toViewport();
		glDisable(GL_BLEND);
	}
	

	shader->disable();


	//irradiance
	//applyIrradiance();


	glEnable(GL_DEPTH_TEST);


	if (show_probes) {
		for (int i = 0; i < probes.size(); ++i)			
			renderProbe(probes[i]);
	}

	light_fbo->unbind();


	if (show_gbuffers)
	{
		glViewport(0, size.y / 2, size.x / 2, size.y / 2);
		gbuffers_fbo->color_textures[0]->toViewport();
		glViewport(size.x / 2, size.y / 2, size.x / 2, size.y / 2);
		gbuffers_fbo->color_textures[1]->toViewport();
		glViewport(0,0, size.x / 2, size.y / 2);
		gbuffers_fbo->color_textures[2]->toViewport();
		glViewport(size.x / 2,0, size.x / 2, size.y / 2);
		shader = GFX::Shader::getDefaultShader("linear_depth");
		shader->enable();
		shader->setUniform("u_camera_nearfar", vec2(camera->near_plane, camera->far_plane));
		gbuffers_fbo->depth_texture->toViewport(shader);
		glViewport(0,0, size.x, size.y);
		shader->disable();
	}
	else
	{
		shader = GFX::Shader::Get("tonemapper");
		shader->enable();
		shader->setUniform("u_scale", tonemapper_scale); //color scale before tonemapper
		shader->setUniform("u_average_lum",average_lum);
		shader->setUniform("u_lumwhite2", lumwhite2);
		shader->setUniform("u_igamma",1.0f/igamma); //inverse gamma

		light_fbo->color_textures[0]->toViewport(shader);

	}
	if (show_volumetric)
		volumetric_fbo->color_textures[0]->toViewport();

}

void SCN::Renderer::captureIrradiance()
{


	//when computing the probes position…

	//define the corners of the axis aligned grid
	//this can be done using the boundings of our scene
	vec3 start_pos(-300, 5, -400);
	vec3 end_pos(300, 150, 400);

	//define how many probes you want per dimension
	vec3 dim(10, 4, 10);

	//compute the vector from one corner to the other
	vec3 delta = (end_pos - start_pos);

	//and scale it down according to the subdivisions
	//we substract one to be sure the last probe is at end pos
	delta.x /= (dim.x - 1);
	delta.y /= (dim.y - 1);
	delta.z /= (dim.z - 1);

	//now delta give us the distance between probes in every axis
	probes.resize(dim.x * dim.y * dim.z);
	//lets compute the centers
	//pay attention at the order at which we add them
	for (int z = 0; z < dim.z; ++z)
		for (int y = 0; y < dim.y; ++y)
			for (int x = 0; x < dim.x; ++x)
			{
				sProbe p;
				p.local.set(x, y, z);

				//index in the linear array
				p.index = x + y * dim.x + z * dim.x * dim.y;

				//and its position
				p.pos = start_pos + delta * vec3(x, y, z);
				probes[p.index] = p;
			}


	//now compute the coeffs for every probe
	for (int iP = 0; iP < probes.size(); ++iP)
	{
		int probe_index = iP;
		sProbe& p = probes[iP];
		captureProbe(p);
		//...
	}


	FILE* f = fopen("irradiance_cache.bin", "wb");
	if (f == NULL)
		return;

	irrCacheInfo.dims = dim;
	irrCacheInfo.start = start_pos;
	irrCacheInfo.end = end_pos;
	irrCacheInfo.num_probes = probes.size();
	fwrite(&irrCacheInfo, sizeof(irrCacheInfo), 1, f);
	fwrite(&probes[0], sizeof(sProbe), probes.size(), f);

	fclose(f);
	uploadIrradianceCache();
}

void SCN::Renderer::loadIrradianceCache()
{
	FILE* f = fopen("irradiance_cache.bin", "rb");
	if (f == NULL)
		return;

	fread(&irrCacheInfo, sizeof(sIrradianceCacheInfo), 1, f);
	probes.resize(irrCacheInfo.num_probes);
	fread(&probes[0], sizeof(sProbe), irrCacheInfo.num_probes, f);

	fclose(f);

	uploadIrradianceCache();
}


void SCN::Renderer::captureProbe(sProbe& probe)
{
	render_mode = eRenderMode::LIGHTMULTIPASS;

	FloatImage images[6]; //here we will store the six views
	Camera cam;
	//set the fov to 90 and the aspect to 1
	cam.setPerspective(90, 1, 0.1, 1000);

	if (!irr_fbo)
	{
		irr_fbo = new GFX::FBO();
		irr_fbo->create(64, 64, 1, GL_RGB, GL_FLOAT, false);
	}

	for (int i = 0; i < 6; ++i) //for every cubemap face
	{
		//compute camera orientation using defined vectors
		vec3 eye = probe.pos;
		vec3 front = cubemapFaceNormals[i][2];
		vec3 center = probe.pos + front;
		vec3 up = cubemapFaceNormals[i][1];
		cam.lookAt(eye, center, up);
		cam.enable();

		//render the scene from this point of view
		irr_fbo->bind();
		renderForward(scene, &cam);
		irr_fbo->unbind();

		//read the pixels back and store in a FloatImage
		images[i].fromTexture(irr_fbo->color_textures[0]);
	}
	render_mode = eRenderMode::DEFERRED;
	//compute the coefficients given the six images
	probe.sh = computeSH(images);

}

void SCN::Renderer::renderProbe(sProbe& probe)
{
	Camera* camera = Camera::current;
	GFX::Shader* shader = GFX::Shader::Get("spherical_probe");
	shader->enable();

	glEnable(GL_CULL_FACE);
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	Matrix44 model;
	model.setTranslation(probe.pos.x, probe.pos.y, probe.pos.z);
	model.scale(10, 10, 10);

	shader->setUniform("u_model", model);
	shader->setUniform3Array("u_coeffs", probe.sh.coeffs[0].v, 9);

	cameraToShader(camera, shader);


	sphere.render(GL_TRIANGLES);
}
void SCN::Renderer::uploadIrradianceCache()
{
	if (probes_texture)
		delete probes_texture;

	vec3 dim = irrCacheInfo.dims;

	//create the texture to store the probes (do this ONCE!!!)
	probes_texture = new GFX::Texture(
		9, //9 coefficients per probe
		probes.size(), //as many rows as probes
		GL_RGB, //3 channels per coefficient
		GL_FLOAT); //they require a high range

		//we must create the color information for the texture. because every SH are 27 floats in the RGB,RGB,... order, we can create an array of SphericalHarmonics and use it as pixels of the texture
	SphericalHarmonics* sh_data = NULL;
	int sizeofdata = dim.x * dim.y * dim.z;
	sh_data = new SphericalHarmonics[sizeofdata];

	//here we fill the data of the array with our probes in x,y,z order
	for (int i = 0; i < probes.size(); ++i)
		sh_data[i] = probes[i].sh;

	//now upload the data to the GPU as a texture
	probes_texture->upload(GL_RGB, GL_FLOAT, false, (uint8*)sh_data);

	//disable any texture filtering when reading
	probes_texture->bind();
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);

	//always free memory after allocating it!!!
	delete[] sh_data;

}

void SCN::Renderer::applyIrradiance()
{
	if (!probes_texture)
		return;
	GFX::Mesh* mesh = GFX::Mesh::getQuad();
	Camera* camera = Camera::current;
	vec2 size = CORE::getWindowSize();

	glDisable(GL_DEPTH_TEST);
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE);
	GFX::Shader* shader = GFX::Shader::Get("irradiance");
	shader->enable();

	shader->setTexture("u_normal_texture", gbuffers_fbo->color_textures[1], 1);
	shader->setTexture("u_depth_texture", gbuffers_fbo->depth_texture, 3);
	shader->setTexture("u_albedo_texture", gbuffers_fbo->color_textures[0], 0);
	shader->setUniform("u_iRes", vec2(1 / size.x, 1 / size.y));
	shader->setUniform("u_ivp", camera->inverse_viewprojection_matrix);
	shader->setUniform("u_camera_position", camera->eye);

	//compute the vector from one corner to the other
	vec3 delta = (irrCacheInfo.end - irrCacheInfo.start);

	//and scale it down according to the subdivisions
	//we substract one to be sure the last probe is at end pos
	delta.x /= (irrCacheInfo.dims.x - 1);
	delta.y /= (irrCacheInfo.dims.y - 1);
	delta.z /= (irrCacheInfo.dims.z - 1);

	//probes
	shader->setUniform("u_irr_start", irrCacheInfo.start);
	shader->setUniform("u_irr_end", irrCacheInfo.end);
	shader->setUniform("u_irr_dims", irrCacheInfo.dims);
	shader->setUniform("u_irr_delta", delta);

	shader->setUniform("u_irr_normal_distance", .1f);
	shader->setTexture("u_probe_texture", probes_texture, 4);
	shader->setUniform("u_num_probes", irrCacheInfo.num_probes);
	shader->setUniform("u_num_probes", irrCacheInfo.num_probes);
	shader->setUniform("u_irr_multiplier", irradiance_multiplier);

	mesh->render(GL_TRIANGLES);


}

void Renderer::debug_shadowmap(){
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_BLEND);

	int x = 310;

	for (auto light : lights) {
		if (!light->shadowmap)
			continue;

		GFX::Shader* shader = GFX::Shader::getDefaultShader("linear_depth");
		shader->enable();
		shader->setUniform("u_camera_nearfar",vec2(light->near_distance,light->max_distance));
		glViewport(x, 100, 256, 256);
		light->shadowmap->toViewport(shader);
		x += 256;
	}

	vec2 size = CORE::getWindowSize();
	glViewport(0, 0, size.x, size.y);
}

void Renderer::generateShadowMaps() {

	Camera camera;
	eRenderMode prev = render_mode;
	render_mode = eRenderMode::FLAT;

	for (auto light : lights) {
		if (!light->cast_shadows)
			continue;
		if (light->light_type == eLightType::POINT || light->light_type == eLightType::NO_LIGHT)
			continue;
		//check if light inside camera
		//TODO

		if (!light->shadowmap_fbo) {
			light->shadowmap_fbo = new GFX::FBO();
			light->shadowmap_fbo->setDepthOnly(1024, 1024);
			light->shadowmap = light->shadowmap_fbo->depth_texture;
		}
		vec3 pos = light->root.model.getTranslation();
		vec3 front = light->root.model.rotateVector(vec3(0,0,-1));
		vec3 up = vec3(0, 1, 0);
		camera.lookAt(pos, pos+front, up);
		if (light->light_type == eLightType::DIRECTIONAL) {
			float halfarea = light->area / 2;
			camera.setOrthographic(-halfarea, halfarea,
				halfarea , -halfarea,
				0.1, light->max_distance);

		}
		else if(light->light_type == eLightType::SPOT) {
			camera.setPerspective(light->cone_info.y * 2, 1.0, light->near_distance, light->max_distance);
		}



		//if (prev == eRenderMode::LIGHTMULTIPASS) 
		//{
		light->shadowmap_fbo->bind();
		renderFrame(scene,&camera);
		light->shadowmap_fbo->unbind();

		//}
		/*
		else if (prev == eRenderMode::LIGHTSINGLEPASS) 
		
		vec4 shadow_region(1024*floor(i/ max_shadowmaps_pow2), 1024*(i% max_shadowmaps_pow2), 1024, 1024);
		
		glBindFramebuffer(GL_READ_FRAMEBUFFER, light->shadowmap_fbo->fbo_id);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, shadow_atlas_fbo->fbo_id);

		// Set the source and destination rectangles for the blit
		glBlitFramebuffer(0, 0, 1024, 1024, // source rectangle
			shadow_region.x, shadow_region.y, shadow_region.x + 1024, shadow_region.y + 1024, // destination rectangle
			GL_DEPTH_BUFFER_BIT, GL_NEAREST);

		// Unbind the FBOs
		glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
		glBindFramebuffer(GL_DRAW_FRAMEBUFFER, 0);
/*
		}
		*/

		light->shadow_viewproj = camera.viewprojection_matrix;


		
	}

	render_mode = prev;

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
	ImGui::Checkbox("Show ShadowMap", &show_shadowmap);
	ImGui::Checkbox("Show shadows", &show_shadows);
	ImGui::Checkbox("Normal maps", &normal_maps);
	//ImGui::SliderInt("Lights shininess coeficient", &shininess_coef, 10, 100);
	ImGui::SliderFloat("Tonemapper scale coeficient", &tonemapper_scale, 0, 2);
	ImGui::SliderFloat("Agerage lumminosity coeficient", &average_lum, 0, 2);
	ImGui::SliderFloat("Lumwhite2 coeficient", &lumwhite2, 0, 2);
	ImGui::SliderFloat("Igamma coeficient", &igamma, 0, 2);

	ImGui::Combo("Render Mode", (int*) & render_mode, "FLAT\0TEXTURED\0LIGHTSMULTIPASS\0LIGHTSSINGLEPASS\0DEFERRED\0", 2);
	ImGui::Checkbox("Show volumetric", &show_volumetric);
	ImGui::DragFloat("Air density", &air_density, 0.0001, 0.0001, 0.1);
	if(render_mode==eRenderMode::DEFERRED)
		ImGui::Checkbox("Show gbuffers", &show_gbuffers);
		ImGui::Checkbox("Show probes", &show_probes);

	if (ImGui::Button("Update Probes"))
	{
		captureIrradiance();
	}
	ImGui::SameLine();
	if (ImGui::Button("Load Probes"))
	{
		loadIrradianceCache();
	}

}

#else
void Renderer::showUI() {}
#endif