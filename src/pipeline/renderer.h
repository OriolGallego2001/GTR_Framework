#pragma once
#include "scene.h"
#include "prefab.h"
#include "../gfx/sphericalharmonics.h"

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
		TEXTURED,
		LIGHTMULTIPASS,
		LIGHTSINGLEPASS,
		DEFERRED
	};

	// This class is in charge of rendering anything in our system.
	// Separating the render from anything else makes the code cleaner
	class Renderer
	{
	public:
		bool render_wireframe;
		bool render_boundaries;
		eRenderMode render_mode;

		//Debug parameters
		bool show_shadowmap;

		GFX::Texture* skybox_cubemap;

		SCN::Scene* scene;

		std::vector<LightEntity*> lights;
		std::vector<LightEntity*> visible_lights;


		//updated every frame
		Renderer(const char* shaders_atlas_filename );

		//just to be sure we have everything ready for the rendering
		void setupScene(Camera* camera);

		//add here your functions
		//...
		class RenderCall {
			public:
				GFX::Mesh* mesh;
				Material* material;
				Matrix44 model;
				
				float distance_to_camera;
		};

		std::vector<RenderCall> opaqueRenderCalls;
		std::vector<RenderCall> transparentRenderCalls;
		std::vector<DecalEntity*> decals;


		bool show_shadows;
		bool normal_maps;
		int shininess_coef;
		bool show_gbuffers;
		bool show_probes;

		GFX::FBO* shadow_atlas_fbo;
		int max_shadowmaps_pow2 = 2;

		GFX::FBO* light_fbo;
		GFX::FBO* ssao_fbo;
		GFX::FBO* irr_fbo;
		GFX::FBO* volumetric_fbo;
		GFX::FBO* gbuffers_fbo;

		float air_density;
		bool show_volumetric;

		struct sProbe {
			vec3 pos; //where is located
			vec3 local; //its ijk pos in the matrix
			int index; //its index in the linear array
			SphericalHarmonics sh; //coeffs
		};
		std::vector<sProbe> probes;

		
		struct sIrradianceCacheInfo {
			int num_probes;
			vec3 dims;
			vec3 start;
			vec3 end;
		};
		sIrradianceCacheInfo irrCacheInfo;
		float irradiance_multiplier;


		float tonemapper_scale; //color scale before tonemapper
		float average_lum;
		float lumwhite2;
		float igamma; //inverse gamma

		GFX::Texture* current_shadow_texture;
		GFX::Texture* cloned_depth_texture;

		float distance(vec3 p1, vec3 p2);
		
		


		//renders several elements of the scene
		void renderScene(SCN::Scene* scene, Camera* camera);
		void renderFrame(SCN::Scene* scene, Camera* camera);
		void renderForward(SCN::Scene* scene, Camera* camera);
		void renderDeferred(Scene* scene, Camera* camera);
		void renderRenderCalls(SCN::Scene* scene, Camera* camera);
		void renderDeferredVolumetric(Scene* scene, Camera* camera);
		void renderWithDecals(Scene* scene, Camera* camera);

		void lightToShader(LightEntity* light, GFX::Shader* shader);



		void renderMeshWithMaterialGBuffers(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void  generateShadowMaps();

		//render the skybox
		void renderSkybox(GFX::Texture* cubemap, float skybox_instensity);
	
		//to render one node from the prefab and its children
		void renderNode(SCN::Node* node, Camera* camera);

		//to render one mesh given its material and transformation matrix
		void renderMeshWithMaterial(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);
		//Add all the lights of the scene
		void renderMeshWithMaterialLight(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void renderMeshWithMaterialFlat(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void renderMeshWithMaterialLightSinglePass(const Matrix44 model, GFX::Mesh* mesh, SCN::Material* material);

		void captureProbe(sProbe& probe);
		void renderProbe(sProbe& probe);
		void captureIrradiance();
		void loadIrradianceCache();
		void applyIrradiance();
		void uploadIrradianceCache();
		GFX::Texture* probes_texture;


		void showUI();

		void cameraToShader(Camera* camera, GFX::Shader* shader); //sends camera uniforms to shader

		//Debug
		void debug_shadowmap();
	};

};