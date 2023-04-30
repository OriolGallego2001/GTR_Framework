#pragma once
#include "scene.h"
#include "prefab.h"

#include "light.h"


//forward declarations
class Camera;
class Skeleton;
namespace GFX {
	class Shader;
	class Mesh;
	class FBO;
}

namespace SCN {

	class Prefab;
	class Material;

	enum eRenderMode {
		FLAT,
		LIGHTMULTIPASS,
		LIGHTSINGLEPASS
	};

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
	public:
		bool render_wireframe;
		bool render_boundaries;
		eRenderMode render_mode;

		GFX::Texture* skybox_cubemap;

		SCN::Scene* scene;

		std::vector<LightEntity*> lights;
		std::vector<LightEntity*> visible_lights;


		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene();

		//add here your functions
		//...
		class RenderCall {
			public:
				GFX::Mesh* mesh;
				Material* material;
				Matrix44 model;
				
				float distance_to_camera;
		};

		std::vector<RenderCall> renderCalls;

		float distance(vec3 p1, vec3 p2);
		


		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);
		//Add all the lights of the scene
		void renderMeshWithMaterialLight(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void renderMeshWithMaterialLightSinglePass(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void showUI();

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader
	};

};