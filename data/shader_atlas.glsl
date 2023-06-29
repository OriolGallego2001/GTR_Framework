//example of some shaders compiled
flat basic.vs flat.fs
texture basic.vs texture.fs
light basic.vs light.fs
lightSinglePass basic.vs lightSinglePass.fs
skybox basic.vs skybox.fs
depth quad.vs depth.fs
multi basic.vs multi.fs
gbuffers basic.vs gbuffers.fs
tonemapper quad.vs tonemapper.fs
ssao quad.vs ssao.fs
decal basic.vs decal.fs
volumetric quad.vs volumetric.fs

deferred_global quad.vs deferred_global.fs
deferred_light basic.vs deferred_light.fs

spherical_probe basic.vs spherical_probe.fs
irradiance quad.vs irradiance.fs
fx_color_correction quad.vs color_correction.fs
fx_blur quad.vs blur.fs
depthOfField quad.vs depthOfField.fs

\basic.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;
in vec4 a_color;

uniform vec3 u_camera_pos;

uniform mat4 u_model;
uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;
out vec4 v_color;

uniform float u_time;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( v_position, 1.0) ).xyz;
	
	//store the color in the varying var to use it from the pixel shader
	v_color = a_color;

	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}

\quad.vs

#version 330 core

in vec3 a_vertex;
in vec2 a_coord;
out vec2 v_uv;

void main()
{	
	v_uv = a_coord;
	gl_Position = vec4( a_vertex, 1.0 );
}


\flat.fs

#version 330 core

uniform vec4 u_color;

out vec4 FragColor;

void main()
{
	FragColor = u_color;
}


\texture.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, v_uv );

	if(color.a < u_alpha_cutoff)
		discard;

	FragColor = color;
}

\tonemapper.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_texture;

out vec4 FragColor;

uniform float u_scale; //color scale before tonemapper
uniform float u_average_lum;
uniform float u_lumwhite2;
uniform float u_igamma; //inverse gamma


void main() {
	vec4 color = texture2D(u_texture, v_uv);
	vec3 rgb = color.xyz;

	float lum = dot(rgb, vec3(0.2126, 0.7152, 0.0722));
	float L = (u_scale / u_average_lum) * lum;
	float Ld = (L * (1.0 + L / u_lumwhite2)) / (1.0 + L);

	rgb = (rgb / lum) * Ld;
	rgb = max(rgb, vec3(0.001));
	rgb = pow(rgb, vec3(u_igamma));
	FragColor = vec4(rgb, color.a);
}

\gbuffers.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

uniform vec4 u_color;
uniform sampler2D u_albedo_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

uniform sampler2D u_emissive_texture;
uniform vec3 u_emissive_factor;
uniform sampler2D u_metalness_texture;

uniform float u_roughness;
uniform float u_metalness;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;
layout(location = 2) out vec4 ExtraColor;

vec3 degamma(vec3 c)
{
	return pow(c, vec3(2.2));
}

vec3 gamma(vec3 c)
{
	return pow(c, vec3(1.0 / 2.2));
}


void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	vec4 alb_tex = texture(u_albedo_texture, v_uv);
	color *= alb_tex;

	if (color.a < u_alpha_cutoff)
		discard;

	vec3 emissive = u_emissive_factor * texture(u_emissive_texture, v_uv).xyz;
	vec2 metalness_roughness = texture(u_metalness_texture, v_uv).yz;


	FragColor = vec4(color.xyz, 1.0);


	NormalColor = vec4(normalize(v_normal)*0.5+vec3(0.5), u_metalness* metalness_roughness.x);
	ExtraColor = vec4(emissive, u_roughness* metalness_roughness.y);
}

\ssao.fs

#version 330 core

in vec2 v_uv;

uniform vec2 u_iRes;
uniform mat4 u_inv_viewproj;
uniform sampler2D u_depth_texture;
uniform sampler2D u_normal_texture;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;
layout(location = 2) out vec4 ExtraColor;

void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;


	float depth = texture(u_depth_texture, uv).x;

	if (depth == 1)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_inv_viewproj * screen_pos;
	vec3 world_pos = proj_worldpos.xyz / proj_worldpos.w;



}


\light.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

//Material prop
uniform vec4 u_color;
uniform sampler2D u_albedo_texture;
uniform sampler2D u_emissive_texture;
uniform vec3 u_emissive_factor;
uniform vec3 u_ambient_light;

uniform vec3 u_camera_position;
uniform float u_roughness;
uniform float u_metalness;


uniform sampler2D u_normalmap_texture;
//uniform vec2 u_uv_scale;

uniform sampler2D u_occlusion_texture;


//Global prop
uniform float u_time;
uniform float u_alpha_cutoff;

#include "lights"


mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx(p);
	vec3 dp2 = dFdy(p);
	vec2 duv1 = dFdx(uv);
	vec2 duv2 = dFdy(uv);

	// solve the linear system
	vec3 dp2perp = cross(dp2, N);
	vec3 dp1perp = cross(N, dp1);
	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

	// construct a scale-invariant frame 
	float invmax = inversesqrt(max(dot(T, T), dot(B, B)));
	return mat3(T * invmax, B * invmax, N);
}

// assume N, the interpolated vertex normal and 
// WP the world position
//vec3 normal_pixel = texture2D( normalmap, uv ).xyz; 
vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255. / 127. - 128. / 127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color_degamma = u_color;
	vec4 albedo = color_degamma;
	vec4 alb_tex = texture( u_albedo_texture, v_uv );
	albedo *= alb_tex;
	

	vec4 occlusion = texture(u_occlusion_texture, v_uv);
	float occlusion_factor = occlusion.r;

	float metalness = occlusion.y * u_metalness;
	float roughness = occlusion.z * u_roughness;
	if(albedo.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	vec3 normal_pixel = texture(u_normalmap_texture, v_uv).rgb;
	vec3 new_normal = perturbNormal(N, v_world_position, v_uv, normal_pixel);
	N = new_normal;


	vec3 light = vec3(0.0);
	float shadow_factor = 1.0;

	if (u_shadow_info.x != 0.0)
		shadow_factor = testShadow(v_world_position);

	if (int(u_light_info.x) == POINTLIGHT || int(u_light_info.x) == SPOTLIGHT)
	{
		vec3 L = u_light_position - v_world_position;
		float dist = length(L);
		L /= dist;
		L = normalize(L);

		float NdotL =  dot(N, L);
		float att = max((u_light_info.z - dist) / u_light_info.z, 0.0);


		vec3 V = normalize(u_camera_position - v_world_position);
		vec3 H = normalize(L + V);

		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix(vec3(0.5), color_degamma.xyz, metalness);
		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness) * color_degamma.xyz;

		float NoV = dot(N, V);
		float LoH = dot(L, H);

		vec3 Fr_d = specularBRDF(roughness, f0, dot(N, H), NoV, dot(N, L), LoH);
		vec3 Fd_d = max(NdotL, 0.0) * diffuseColor;
		vec3 direct = Fr_d + Fd_d;


		if (int(u_light_info.x) == SPOTLIGHT) {
			float cos_angle = dot(u_light_front, L);
			if (cos_angle < u_light_cone.y)
				att = 0;
			else if (cos_angle < u_light_cone.x)
				att *= 1 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
		}


		vec3 lightParams = u_light_color * att * shadow_factor;
		light += lightParams * direct;
	}
	else if (int(u_light_info.x) == DIRECTIONAL) {
		float NdotL = dot(N, u_light_front);
		//light += max(NdotL, 0.0) * u_light_color * shadow_factor;

		vec3 L = u_light_front;
		vec3 V = normalize(u_camera_position - v_world_position);
		vec3 H = normalize(L + V);

		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix(vec3(0.5), color_degamma.xyz, metalness);
		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness) * color_degamma.xyz;

		float NoV = max(dot(N, V),0);
		float LoH = max(dot(L, H),0);

		vec3 Fr_d = specularBRDF(roughness, f0, max(dot(N, H),0), NoV, max(dot(N, L),0), LoH);
		vec3 Fd_d = max(NdotL, 0.0) * diffuseColor;
		vec3 direct = Fr_d + Fd_d;

		vec3 lightParams = u_light_color * shadow_factor;
		light += lightParams * direct;
	}

	light += u_ambient_light * occlusion_factor;

	vec3 color = albedo.xyz * light;

	color += u_emissive_factor * texture(u_emissive_texture, v_uv).xyz;


	FragColor = vec4(color,albedo.a);
}

\lights

#define NOLIGHT 0
#define POINTLIGHT 1
#define SPOTLIGHT 2
#define DIRECTIONAL 3

uniform vec3 u_light_position;
uniform vec3 u_light_color;
uniform vec4 u_light_info; //(light_type, near_distance, max_distance, light shininess)
uniform vec3 u_light_front;
uniform vec2 u_light_cone; //(cos(min_angle), cos(max_angle))


uniform vec2 u_shadow_info; // 0 or 1 if it has shadowmap or not; bias
uniform sampler2D u_shadowmap;
uniform mat4 u_shadow_viewproj;

#define PI 3.14159265359

vec3 degamma(vec3 c)
{
	return pow(c, vec3(2.2));
}

vec3 gamma(vec3 c)
{
	return pow(c, vec3(1.0 / 2.2));
}


float D_GGX(const in float NoH,
	const in float linearRoughness)
{
	float a2 = linearRoughness * linearRoughness;
	float f = (NoH * NoH) * (a2 - 1.0) + 1.0;
	return a2 / (PI * f * f);
}

float GGX(float NdotV, float k) {
	return NdotV / (NdotV * (1.0 - k) + k);
}

float G_Smith(float NdotV, float NdotL, float roughness)
{
	float k = pow(roughness + 1.0, 2.0) / 8.0;
	return GGX(NdotL, k) * GGX(NdotV, k);
}


// Fresnel term with colorized fresnel
vec3 F_Schlick( const in float VoH, 
const in vec3 f0)
{
	float f = pow(1.0 - VoH, 5.0);
	return f0 + (vec3(1.0) - f0) * f;
}


vec3 specularBRDF(float roughness, vec3 f0,
	float NoH, float NoV, float NoL, float LoH)
{
	float a = roughness * roughness;

	// Normal Distribution Function
	float D = D_GGX(NoH, a);

	// Fresnel Function
	vec3 F = F_Schlick(LoH, f0);

	// Visibility Function (shadowing/masking)
	float G = G_Smith(NoV, NoL, roughness);

	// Norm factor
	vec3 spec = D * G * F;
	spec /= (4.0 * NoL * NoV + 1e-6);

	return spec;
}


float testShadow(vec3 pos)
{
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj * vec4(pos, 1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = proj_pos.xy / proj_pos.w;

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);

	//it is outside on the sides
	if (shadow_uv.x < 0.0 || shadow_uv.x > 1.0 ||
		shadow_uv.y < 0.0 || shadow_uv.y > 1.0)
		return 0.0;

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_shadow_info.y) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;


	//read depth from depth buffer in [0..+1] non-linear
	float shadow_depth = texture(u_shadowmap, shadow_uv).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if (shadow_depth < real_depth)
		shadow_factor = 0.0;
	return shadow_factor;
}


\lightSinglePass.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;
in vec4 v_color;

//Material prop
uniform vec4 u_color;
uniform sampler2D u_albedo_texture;
uniform sampler2D u_emissive_texture;
uniform vec3 u_emissive_factor;
uniform sampler2D u_occlusion_texture;
uniform sampler2D u_normalmap_texture;
uniform vec3 u_camera_position;

//Global prop
uniform float u_time;
uniform float u_alpha_cutoff;

uniform vec3 u_ambient_light;

#define NOLIGHT 0
#define POINTLIGHT 1
#define SPOTLIGHT 2
#define DIRECTIONAL 3

const int MAX_LIGHTS = 4;
uniform vec4 u_light_info[MAX_LIGHTS]; //(light_type, near_distance, max_distance, light_shininess)
uniform vec3 u_light_position[MAX_LIGHTS];
uniform vec3 u_light_color[MAX_LIGHTS];
uniform vec3 u_light_front[MAX_LIGHTS];
uniform vec2 u_light_cone[MAX_LIGHTS];
uniform int u_num_lights;

uniform sampler2D u_shadow_atlas;
uniform vec3 u_light_atlas_frame[MAX_LIGHTS]; //[framex0, framey0, bias]
uniform mat4 u_shadow_viewproj[MAX_LIGHTS];

out vec4 FragColor;

mat3 cotangent_frame(vec3 N, vec3 p, vec2 uv)
{
	// get edge vectors of the pixel triangle
	vec3 dp1 = dFdx(p);
	vec3 dp2 = dFdy(p);
	vec2 duv1 = dFdx(uv);
	vec2 duv2 = dFdy(uv);

	// solve the linear system
	vec3 dp2perp = cross(dp2, N);
	vec3 dp1perp = cross(N, dp1);
	vec3 T = dp2perp * duv1.x + dp1perp * duv2.x;
	vec3 B = dp2perp * duv1.y + dp1perp * duv2.y;

	// construct a scale-invariant frame 
	float invmax = inversesqrt(max(dot(T, T), dot(B, B)));
	return mat3(T * invmax, B * invmax, N);
}

// assume N, the interpolated vertex normal and 
// WP the world position
//vec3 normal_pixel = texture2D( normalmap, uv ).xyz; 
vec3 perturbNormal(vec3 N, vec3 WP, vec2 uv, vec3 normal_pixel)
{
	normal_pixel = normal_pixel * 255. / 127. - 128. / 127.;
	mat3 TBN = cotangent_frame(N, WP, uv);
	return normalize(TBN * normal_pixel);
}



float testShadow(vec3 pos, int i)
{
	//project our 3D position to the shadowmap
	vec4 proj_pos = u_shadow_viewproj[i] * vec4(pos, 1.0);

	//from homogeneus space to clip space
	vec2 shadow_uv = proj_pos.xy / proj_pos.w;

	//from clip space to uv space
	shadow_uv = shadow_uv * 0.5 + vec2(0.5);

	//it is outside on the sides
	if (shadow_uv.x < 0.0 || shadow_uv.x > 1.0 ||
		shadow_uv.y < 0.0 || shadow_uv.y > 1.0)
		return 0.0;

	//get point depth [-1 .. +1] in non-linear space
	float real_depth = (proj_pos.z - u_light_atlas_frame[i].z) / proj_pos.w;

	//normalize from [-1..+1] to [0..+1] still non-linear
	real_depth = real_depth * 0.5 + 0.5;

	shadow_uv = shadow_uv + u_light_atlas_frame[i].xy;
	//read depth from depth buffer in [0..+1] non-linear
	float shadow_depth = texture(u_shadow_atlas, shadow_uv).x;

	//compute final shadow factor by comparing
	float shadow_factor = 1.0;

	//we can compare them, even if they are not linear
	if (shadow_depth < real_depth)
		shadow_factor = 0.0;
	return shadow_factor;
}

void main()
{
	vec2 uv = v_uv;
	vec4 albedo = u_color;
	albedo *= texture( u_albedo_texture, v_uv );
	if(albedo.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);
	vec3 light = vec3(0.0);

	vec3 normal_pixel = texture(u_normalmap_texture, v_uv).rgb;
	vec3 new_normal = perturbNormal(N, v_world_position, v_uv, normal_pixel);
	N = new_normal;

	vec4 occlusion = texture(u_occlusion_texture, v_uv);
	float occlusion_factor = occlusion.r;


	int atlas_frame = 0;
	for( int i = 0; i < u_num_lights; ++i )
	{
		float shadow_factor = 1.0;
			
		if(int(u_light_info[i].x) == POINTLIGHT || int(u_light_info[i].x) == SPOTLIGHT)
		{
			vec3 L =  u_light_position[i] - v_world_position;
			float dist = length(L);
			L /= dist;

			float NdotL = dot(N, L);
			float att =max( (u_light_info[i].z-dist)/u_light_info[i].z, 0.0);

			// Compute specular light
			vec3 V = normalize(u_camera_position - v_world_position);
			vec3 R = reflect(N, V);
			float specular = pow(max(dot(R, V), 0.0), u_light_info[i].w);
			vec3 specular_color = u_light_color[i] * u_color.xyz * u_color.w;
			light += specular * specular_color;


				
			if( int(u_light_info[i].x) == SPOTLIGHT){

				shadow_factor = testShadow(v_world_position,i);
				atlas_frame++;

				float cos_angle = dot(u_light_front[i], L);
				if(cos_angle<u_light_cone[i].y)
					att=0;
				else if(cos_angle<u_light_cone[i].x)
					att *= 1-  (cos_angle - u_light_cone[i].x) / (u_light_cone[i].y - u_light_cone[i].x);
			}

			light += max(NdotL, 0.0) * u_light_color[i] * att * shadow_factor;
		}
		else if(int(u_light_info[i].x) == DIRECTIONAL){
			shadow_factor = testShadow(v_world_position,i);
			atlas_frame++;
			float NdotL = dot(N, u_light_front[i]);
			light += max(NdotL, 0.0) * u_light_color[i] * shadow_factor;
		}
		
	}

	light += u_ambient_light * occlusion_factor;
	
	vec3 color = albedo.xyz * light;
	
	color += u_emissive_factor * texture(u_emissive_texture, v_uv).xyz;


	FragColor = vec4(color,albedo.a);
}


\skybox.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;

uniform samplerCube u_texture;
uniform vec3 u_camera_position;
out vec4 FragColor;

void main()
{
	vec3 E = v_world_position - u_camera_position;
	vec4 color = texture( u_texture, E );
	FragColor = color;
}


\multi.fs

#version 330 core

in vec3 v_position;
in vec3 v_world_position;
in vec3 v_normal;
in vec2 v_uv;

uniform vec4 u_color;
uniform sampler2D u_texture;
uniform float u_time;
uniform float u_alpha_cutoff;

layout(location = 0) out vec4 FragColor;
layout(location = 1) out vec4 NormalColor;

void main()
{
	vec2 uv = v_uv;
	vec4 color = u_color;
	color *= texture( u_texture, uv );

	if(color.a < u_alpha_cutoff)
		discard;

	vec3 N = normalize(v_normal);

	FragColor = color;
	NormalColor = vec4(N,1.0);
}


\depth.fs

#version 330 core

uniform vec2 u_camera_nearfar;
uniform sampler2D u_texture; //depth map
in vec2 v_uv;
out vec4 FragColor;

void main()
{
	float n = u_camera_nearfar.x;
	float f = u_camera_nearfar.y;
	float z = texture2D(u_texture,v_uv).x;
	if( n == 0.0 && f == 1.0 )
		FragColor = vec4(z);
	else
		FragColor = vec4( n * (z + 1.0) / (f + n - z * (f - n)) );
}


\instanced.vs

#version 330 core

in vec3 a_vertex;
in vec3 a_normal;
in vec2 a_coord;

in mat4 u_model;

uniform vec3 u_camera_pos;

uniform mat4 u_viewprojection;

//this will store the color for the pixel shader
out vec3 v_position;
out vec3 v_world_position;
out vec3 v_normal;
out vec2 v_uv;

void main()
{	
	//calcule the normal in camera space (the NormalMatrix is like ViewMatrix but without traslation)
	v_normal = (u_model * vec4( a_normal, 0.0) ).xyz;
	
	//calcule the vertex in object space
	v_position = a_vertex;
	v_world_position = (u_model * vec4( a_vertex, 1.0) ).xyz;
	
	//store the texture coordinates
	v_uv = a_coord;

	//calcule the position of the vertex using the matrices
	gl_Position = u_viewprojection * vec4( v_world_position, 1.0 );
}


\deferred_global.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_extra_texture;
uniform sampler2D u_depth_texture;
uniform vec3 u_ambient_light;
#include "lights"

out vec4 FragColor;

void main()
{
	vec2 uv = v_uv;

	float depth = texture(u_depth_texture, uv).x;

	if (depth == 1)
		discard;

	vec4 albedo = texture (u_albedo_texture, uv);

	vec4 extra = texture (u_extra_texture, uv);
	vec4 normal_info = texture (u_normal_texture, uv);

	vec3 ambient = u_ambient_light;
	vec3 N = normalize(normal_info.xyz * 2.0 - vec3(1.0));

	vec4 color = vec4(0.0);
	color.xyz += extra.xyz + ambient * albedo.xyz;

	

	FragColor = color;
	gl_FragDepth = depth;
}

\deferred_light.fs

#version 330 core

in vec2 v_uv;

uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture; // Normal + metalness
uniform sampler2D u_extra_texture; // Emissive + roughness
uniform sampler2D u_depth_texture;
uniform mat4 u_inv_viewproj;
uniform vec2 u_iRes;
uniform mat4 u_model;
uniform vec3 u_camera_position;

uniform vec4 u_color;


#include "lights"

out vec4 FragColor;

void main()
{

	//vec2 uv = v_uv;
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;
	float depth = texture(u_depth_texture, uv).x;

	if (depth == 1)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_inv_viewproj * screen_pos;
	vec3 world_pos = proj_worldpos.xyz / proj_worldpos.w;

	vec4 albedo = texture (u_albedo_texture, uv);

	vec4 extra = texture (u_extra_texture, uv);
	vec4 normal_info = texture (u_normal_texture, uv);
	vec3 N = normalize(normal_info.xyz * 2.0 - vec3(1.0));
	//Linearize colors

	vec3 light_color = u_light_color;

	vec3 light = vec3(0.0);
	float shadow_factor = 1.0;

	if (u_shadow_info.x != 0.0)
		shadow_factor = testShadow(world_pos);

	if (int(u_light_info.x) == DIRECTIONAL)
	{
		float NdotL = dot(N, u_light_front);
		//light += max(NdotL, 0.0) * u_light_color * shadow_factor;
		vec3 L = normalize(u_light_front);

		vec3 V = normalize(u_camera_position - world_pos);
		vec3 H = normalize(L + V);
		float metalness = normal_info.w;
		float roughness = extra.w;

		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix(vec3(0.5), u_color.xyz, metalness);
		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness) * albedo.xyz;

		float NoV = max(dot(N, V), 0);
		float LoH = max(dot(L, H), 0);

		vec3 Fr_d = specularBRDF(roughness, f0, max(dot(N, H), 0), NoV, max(dot(N, L), 0), LoH);

		vec3 Fd_d = diffuseColor * NdotL;

		//add diffuse and specular reflection
		vec3 direct = Fr_d + Fd_d;
		//compute how much light received the pixel
		vec3 lightParams = light_color * shadow_factor;// *max(NdotL, 0.0);

		//modulate direct light by light received
		light += lightParams * direct;
	}
	else if (int(u_light_info.x) == POINTLIGHT || int(u_light_info.x) == SPOTLIGHT)
	{
		vec3 L = u_light_position - world_pos;
		float dist = length(L);
		L /= dist;
		float NdotL = dot(N, L);
		float att = max((u_light_info.z - dist) / u_light_info.z, 0.0);

		vec3 V = normalize(u_camera_position - world_pos);
		vec3 H = normalize(L + V);
		float metalness = normal_info.w;
		float roughness = extra.w;

		//we compute the reflection in base to the color and the metalness
		vec3 f0 = mix(vec3(0.5), u_color.xyz, metalness);
		//metallic materials do not have diffuse
		vec3 diffuseColor = (1.0 - metalness)* albedo.xyz;

		float NoV = max(dot(N, V),0);
		float LoH = max(dot(L, H),0);

		vec3 Fr_d = specularBRDF(roughness, f0, max(dot(N,H),0), NoV, max(dot(N,L),0), LoH);

		vec3 Fd_d = diffuseColor * NdotL;

		//add diffuse and specular reflection
		vec3 direct = Fr_d + Fd_d;

		if (int(u_light_info.x) == SPOTLIGHT) {
			float cos_angle = dot(u_light_front, L);
			if (cos_angle < u_light_cone.y)
				att = 0;
			else if (cos_angle < u_light_cone.x)
				att *= 1 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
		}

		//compute how much light received the pixel
		vec3 lightParams = light_color * att * shadow_factor;// *max(NdotL, 0.0);

		//modulate direct light by light received
		light += lightParams *direct;
	}

	vec4 color = vec4(0.0);
	color.xyz = light * albedo.xyz;

	//color.xyz = mod(abs(world_pos * 0.01),vec3(1.0));



	FragColor = color;
	gl_FragDepth = depth;
}


\decal.fs
#version 330 core

in vec2 v_uv;

uniform mat4 u_inverse_viewprojection;
uniform mat4 u_light_shadowmap_vp;
uniform vec2 u_iRes;
uniform vec3 u_camera_position;

uniform mat4 u_model;
uniform mat4 u_imodel;

uniform sampler2D u_depth_texture;
uniform sampler2D u_decal_texture;
uniform sampler2D u_gb0_texture;
uniform sampler2D u_gb1_texture;
uniform sampler2D u_gb2_texture;




out vec4 FragColor;

void main() {

	vec2 uv = gl_FragCoord.xy * u_iRes.xy;

	float depth = texture(u_depth_texture, uv).x;


	vec4 screen_position = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);

	vec4 proj_worldpos = u_inverse_viewprojection * screen_position;
	vec3 worldPos = proj_worldpos.xyz / proj_worldpos.w;

	vec3 localpos = (u_imodel * vec4(worldPos, 1.0)).xyz;

	if (localpos.x < -1.0 || localpos.x> 1.0 ||
		localpos.y < -1.0 || localpos.y> 1.0 ||
		localpos.z < -1.0 || localpos.z> 1.0)
		discard;

	//vec4 color= vec4(worldPos.xyz*.01,1.0);
	vec2 decal_uv = localpos.xz * 0.5 + vec2(0.5);
	vec4 decalColor = texture(u_decal_texture, decal_uv);
	vec4 color = decalColor;

	FragColor = color;
}


\volumetric.fs

#version 330 core

in vec2 v_uv;

uniform vec2 u_iRes;
uniform mat4 u_ivp;
uniform sampler2D u_depth_texture;
uniform vec3 u_camera_position;


#include "lights"
uniform float u_air_density;

#define SAMPLES 64
out vec4 FragColor;

float rand(vec2 co)
{
	return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453123);
}
vec3 computeVolLight(vec3 world_pos)
{
	vec3 light = vec3(0.0);
	float shadow_factor = testShadow(world_pos);
	if (int(u_light_info.x) == POINTLIGHT || int(u_light_info.x) == SPOTLIGHT)
	{
		vec3 L = u_light_position - world_pos;
		float dist = length(L);
		L /= dist;

		float att = max((u_light_info.z - dist) / u_light_info.z, 0.0);

		if (int(u_light_info.x) == SPOTLIGHT)
		{
			float cos_angle = dot(u_light_front, L);
			if (cos_angle < u_light_cone.y)
				att = 0;
			else if (cos_angle < u_light_cone.x)
				att *= 1 - (cos_angle - u_light_cone.x) / (u_light_cone.y - u_light_cone.x);
		}

		light += u_light_color * att * shadow_factor;
	}
	else if (int(u_light_info.x) == DIRECTIONAL)
	{
		light += u_light_color * shadow_factor;

	}

	return light;
}
void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;

	float depth = texture(u_depth_texture, uv).x;


	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_pos = proj_worldpos.xyz / proj_worldpos.w;
	float dist = length(u_camera_position - world_pos);

	vec3 ray_start = u_camera_position;
	vec3 ray_dir = (world_pos - ray_start);
	float ray_length = length(ray_dir);
	ray_dir /= ray_length;
	ray_dir = normalize(ray_dir);
	ray_length = min(500.0, ray_length); //max ray

	float step_dist = ray_length / float(SAMPLES);

	ray_start += ray_dir * rand(uv);

	vec3 current_pos = ray_start;
	vec3 ray_offset = ray_dir * step_dist;

	vec3 irradiance = vec3(0.0);
	float transparency = 1.0;
	float air_step = u_air_density * step_dist;

	for (int i = 0; i < SAMPLES; ++i)
	{
		//evaluate contribution
		//...

		vec3 light = computeVolLight(current_pos);

		//accumulate the amount of light
		irradiance += light * transparency * air_step;

		irradiance = max(light, irradiance);
		//advance to next position
		current_pos.xyz += ray_offset;
		//reduce visibility
		transparency -= air_step;
		//too dense, nothing can be seen behind
		if (transparency < 0.001)
			break;

	}

	FragColor = vec4(irradiance, 1-transparency);


}


\irradiance.fs

#version 330 core

const float Pi = 3.141592654;
const float CosineA0 = Pi;
const float CosineA1 = (2.0 * Pi) / 3.0;
const float CosineA2 = Pi * 0.25;
struct SH9 { float c[9]; }; //to store weights
struct SH9Color { vec3 c[9]; }; //to store colors

in vec2 v_uv;
uniform sampler2D u_albedo_texture;
uniform sampler2D u_normal_texture;
uniform sampler2D u_depth_texture;
uniform vec3 u_irr_start;
uniform vec3 u_irr_dims;
uniform vec3 u_irr_end;
uniform vec3 u_irr_delta;
uniform float u_irr_normal_distance;
uniform int u_num_probes;
uniform sampler2D u_probes_texture;
uniform mat4 u_ivp;
uniform vec3 u_iRes;
uniform float u_irr_multiplier;

out vec4 FragColor;

void SHCosineLobe(in vec3 dir, out SH9 sh) //SH9
{
	// Band 0
	sh.c[0] = 0.282095 * CosineA0;
	// Band 1
	sh.c[1] = 0.488603 * dir.y * CosineA1;
	sh.c[2] = 0.488603 * dir.z * CosineA1;
	sh.c[3] = 0.488603 * dir.x * CosineA1;
	// Band 2
	sh.c[4] = 1.092548 * dir.x * dir.y * CosineA2;
	sh.c[5] = 1.092548 * dir.y * dir.z * CosineA2;
	sh.c[6] = 0.315392 * (3.0 * dir.z * dir.z - 1.0) * CosineA2;
	sh.c[7] = 1.092548 * dir.x * dir.z * CosineA2;
	sh.c[8] = 0.546274 * (dir.x * dir.x - dir.y * dir.y) * CosineA2;
}

vec3 ComputeSHIrradiance(in vec3 normal, in SH9Color sh)
{
	// Compute the cosine lobe in SH, oriented about the normal direction
	SH9 shCosine;
	SHCosineLobe(normal, shCosine);
	// Compute the SH dot product to get irradiance
	vec3 irradiance = vec3(0.0);
	for (int i = 0; i < 9; ++i)
		irradiance += sh.c[i] * shCosine.c[i];

	return irradiance;
}
vec3 computeIrr(vec3 local_indices, vec3 N)
{
	//compute in which row is the probe stored
	float row = local_indices.x +
		local_indices.y * u_irr_dims.x +
		local_indices.z * u_irr_dims.x * u_irr_dims.y;

	//find the UV.y coord of that row in the probes texture
	float row_uv = (row + 1.0) / (u_num_probes + 1.0);


	SH9Color sh;

	//fill the coefficients
	const float d_uvx = 1.0 / 9.0;
	for (int i = 0; i < 9; ++i)
	{
		vec2 coeffs_uv = vec2((float(i) + 0.5) * d_uvx, row_uv);
		sh.c[i] = texture(u_probes_texture, coeffs_uv).xyz;
	}

	//now we can use the coefficients to compute the irradiance
	return max(ComputeSHIrradiance(N, sh) * u_irr_multiplier, vec3(0.0));
}


void main()
{
	vec2 uv = gl_FragCoord.xy * u_iRes.xy;


	float depth = texture(u_depth_texture, uv).x;

	if (depth == 1.0)
		discard;

	vec4 screen_pos = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_ivp * screen_pos;
	vec3 world_pos = proj_worldpos.xyz / proj_worldpos.w;


	vec4 albedo = texture (u_albedo_texture, uv);
	vec4 normal_info = texture (u_normal_texture, uv);
	vec3 N = normalize(normal_info.xyz * 2.0 - vec3(1.0));

	//computing nearest probe index based on world position
	vec3 irr_range = u_irr_end - u_irr_start;
	vec3 irr_local_pos = clamp(world_pos - u_irr_start
		+ N * u_irr_normal_distance, //offset a little
		vec3(0.0), irr_range);

	//convert from world pos to grid pos
	vec3 irr_norm_pos = irr_local_pos / u_irr_delta;

	//round values as we cannot fetch between rows for now
	vec3 local_indices = floor(irr_norm_pos);

	vec3 factors = irr_norm_pos - local_indices;

	vec3 indicesLBF = local_indices;
	//example for right bottom far index
	vec3 indicesRBF = local_indices;
	indicesRBF.x += 1; //from left to right

	vec3 indicesLTF = local_indices;
	indicesLTF.y += 1;

	vec3 indicesRTF = local_indices;
	indicesRBF.x += 1;
	indicesRBF.y += 1;

	vec3 indicesLBN = local_indices;
	indicesLBN.z += 1;

	vec3 indicesRBN = local_indices;
	indicesRBN.z += 1;
	indicesRBN.x += 1;

	vec3 indicesLTN = local_indices;
	indicesLTN.z += 1;
	indicesLTN.y += 1;

	vec3 indicesRTN = local_indices;
	indicesRTN.x += 1;
	indicesRTN.y += 1;
	indicesRTN.z += 1;


	//compute irradiance for every corner
	vec3 irrLBF = computeIrr(indicesLBF, N);
	vec3 irrRBF = computeIrr(indicesRBF, N);
	vec3 irrLTF = computeIrr(indicesLTF, N);
	vec3 irrRTF = computeIrr(indicesRTF, N);
	vec3 irrLBN = computeIrr(indicesLBN, N);
	vec3 irrRBN = computeIrr(indicesRBN, N);
	vec3 irrLTN = computeIrr(indicesLTN, N);
	vec3 irrRTN = computeIrr(indicesRTN, N);

	vec3 irrTF = mix(irrLTF, irrRTF, factors.x);
	vec3 irrBF = mix(irrLBF, irrRBF, factors.x);
	vec3 irrTN = mix(irrLTN, irrRTN, factors.x);
	vec3 irrBN = mix(irrLBN, irrRBN, factors.x);

	vec3 irrT = mix(irrTF, irrTN, factors.z);
	vec3 irrB = mix(irrBF, irrBN, factors.z);

	vec3 irr = mix(irrB, irrT, factors.y);

	// irr *= albedo.xyz;

	FragColor = vec4( 1.0);


}



\spherical_probe.fs

#version 330 core

in vec3 v_world_position;
in vec3 v_normal;
uniform vec3 u_coeffs[9];


const float Pi = 3.141592654;
const float CosineA0 = Pi;
const float CosineA1 = (2.0 * Pi) / 3.0;
const float CosineA2 = Pi * 0.25;
struct SH9 { float c[9]; }; //to store weights
struct SH9Color { vec3 c[9]; }; //to store colors

void SHCosineLobe(in vec3 dir, out SH9 sh) //SH9
{
	// Band 0
	sh.c[0] = 0.282095 * CosineA0;
	// Band 1
	sh.c[1] = 0.488603 * dir.y * CosineA1;
	sh.c[2] = 0.488603 * dir.z * CosineA1;
	sh.c[3] = 0.488603 * dir.x * CosineA1;
	// Band 2
	sh.c[4] = 1.092548 * dir.x * dir.y * CosineA2;
	sh.c[5] = 1.092548 * dir.y * dir.z * CosineA2;
	sh.c[6] = 0.315392 * (3.0 * dir.z * dir.z - 1.0) * CosineA2;
	sh.c[7] = 1.092548 * dir.x * dir.z * CosineA2;
	sh.c[8] = 0.546274 * (dir.x * dir.x - dir.y * dir.y) * CosineA2;
}

vec3 ComputeSHIrradiance(in vec3 normal, in SH9Color sh)
{
	// Compute the cosine lobe in SH, oriented about the normal direction
	SH9 shCosine;
	SHCosineLobe(normal, shCosine);
	// Compute the SH dot product to get irradiance
	vec3 irradiance = vec3(0.0);
	for (int i = 0; i < 9; ++i)
		irradiance += sh.c[i] * shCosine.c[i];

	return irradiance;
}


out vec4 FragColor;

void main()
{
	vec4 color = vec4(1.0);
	vec3 N = normalize(v_normal);
	SH9Color sh;
	for (int i = 0; i < 9; i++)
		sh.c[i] = u_coeffs[i];
	color.xyz = max(vec3(0.0), ComputeSHIrradiance(N, sh));
	FragColor = color;
}

\color_correction.fs
#version 330 core

uniform sampler2D u_texture;
in vec2 v_uv;
out vec4 FragColor;

uniform float u_brightness;

void main() {
	vec4 color = texture(u_texture, v_uv);
	color *= u_brightness;
	FragColor = color;

}


\blur.fs
#version 330 core

in vec2 v_uv;
out vec4 FragColor;

//linear blur shader of 9 samples
uniform sampler2D u_texture;
uniform vec2 u_offset;
uniform float u_intensity;



void main() {
	vec4 sum = vec4(0.0);
	sum += texture(u_texture, v_uv + u_offset * -4.0) * 0.05 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * -3.0) * 0.09 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * -2.0) * 0.12 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * -1.0) * 0.15 / 0.98;
	sum += texture(u_texture, v_uv) * 0.16 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * 4.0) * 0.05 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * 3.0) * 0.09 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * 2.0) * 0.12 / 0.98;
	sum += texture(u_texture, v_uv + u_offset * 1.0) * 0.15 / 0.98;
	FragColor = u_intensity * sum;
}



\depthOfField.fs
#version 330 core


uniform sampler2D u_depth_texture;
uniform mat4 u_inverse_viewprojection;
uniform sampler2D focusTexture;
uniform sampler2D u_texture;

uniform vec2 focusPoint;
uniform vec2 nearFar;


uniform float u_minDistance;
uniform float u_maxDistance;
uniform float u_alpha_cutoff;

uniform vec2 u_iRes;


out vec4 FragColor;

void main() {
	float far = nearFar.y;

	vec2 uv = gl_FragCoord.xy * u_iRes.xy;

	float depth = texture(u_depth_texture, uv).x;

	vec4 focusColor = texture(focusTexture, uv);
	vec4 outOfFocusColor = texture(u_texture, uv);

	FragColor = focusColor;


	vec4 screen_position = vec4(uv.x * 2.0 - 1.0, uv.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_worldpos = u_inverse_viewprojection * screen_position;
	vec3 worldpos = proj_worldpos.xyz / proj_worldpos.w;

	//if (proj_worldpos.a <= 0) { FragColor = vec4(1.0); return; }
	//if(focusColor.a == 0.0f || outOfFocusColor.a==0.0f )
		  //discard;


	vec4 focusPointPos = vec4(focusPoint.x * 2.0 - 1.0, focusPoint.y * 2.0 - 1.0, depth * 2.0 - 1.0, 1.0);
	vec4 proj_focuspos = u_inverse_viewprojection * focusPointPos;
	vec3 focuspos = proj_focuspos.xyz / proj_focuspos.w;

	float blur =
		smoothstep
		(u_minDistance
			, u_maxDistance
			, length(worldpos - focuspos)
		);

	FragColor = mix(focusColor, outOfFocusColor, blur);



}

