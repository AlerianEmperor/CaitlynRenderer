#ifndef _SCENE_SANG_H_
#define _SCENE_SANG_H_
#include <fstream>
#include "Shader.h"
#include <glm\glm.hpp>
#include "Rnd.h"
#include "Camera.h"
#include "Quad.h"
//#include "bvh.h"
//#include "sbvh2.h"
#include "sbvh.h"
//#include "BVH_Optimize.h"
#include "Geometry.h"
#include "Triangle.h"
#include <unordered_map>
#include <unordered_set>
#include <iterator>//back_insertor
#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"
#define STB_IMAGE_RESIZE_IMPLEMENTATION
#include "stb_image_resize.h"

using namespace std;
using namespace glm;

/*struct Material
{
vec4 color;
vec4 emission;
vec4 param;
vec4 texIDs;
};*/

//int screen_width = 500;
//int screen_height = 500;
//float screen_width = 500, screen_height = 500;
float screen_width = 700, screen_height = 700;

enum RendererType
{
	Renderer_Progressive,
	Renderer_Tiled,
};

struct RenderOptions
{
	RenderOptions()
	{
		rendererType = Renderer_Progressive;
		maxSamples = 1024;
		maxDepth = 3;
		numTilesX = 5;
		numTilesY = 5;
		useEnvMap = false;
		resolution = glm::vec2(700, 700);
		hdrMultiplier = 1.0f;
		require_tex_width = 256;
		require_tex_height = 256;
	}
	//std::string rendererType;
	int rendererType; // see RendererType
	glm::ivec2 resolution;
	int maxSamples;
	int maxDepth;
	int numTilesX;
	int numTilesY;
	bool useEnvMap;
	float hdrMultiplier;
	int require_tex_width;//require_tex_width and require_tex_height was use to transform all texture into the same size
	int require_tex_height;//so they can be store in a 2D texture array and sample in rendering process
};

RenderOptions renderOptions;

struct Material
{
	vec4 albedo;  //xyz, w = alpha value, 0 = complete transparent, 1 = complete opaque
	vec4 emission;//xyz, w = -1 : not emission, w = 0, 1, .... : emission and light_index
	vec4 specular;//xyz, w = 0 : not specular matrial, w = 1 : specular material
	vec4 tex_ind = vec4(-1);//albedo ind, normal_ind, specular_ind, metallic_roughness_ind
	Material()
	{
		tex_ind = vec4(-1);
	}
};



struct Texture
{
	vector<vec3> albedo;
	vector<vec3> normal;
	vector<vec3> metallic_roughness;

	//store every texture in one place might speed thing up but can cause overflow
	//unsigned char* albedo;
	//unsigned char* normal;
	//unsigned char* metallic_roughness;


	int num_albedo = 0;
	int num_normal = 0;
	int num_metallic = 0;

	vec2 albedo_size = vec2(512, 512);
	vec2 normal_size;
	vec2 metallic_size;
	Texture() {}
};

enum MaterialType
{
	Diffuse_type,
	Mirror_type,
	Glass_type,
	Glass_Color_type,
	Glass_No_Refract_type,
	Rough_Dielectric_type,
	Conductor_type,
	RoughConductor_type,
	RoughConductorComplex_type,
	RoughConductorSimple_type,
	Plastic_type,
	RoughPlastic_type,
	RoughPlastic_Specular_type,
	ThinSheet_type,
	ThinDielectric_type,
	SmoothCoat_type,
	//Transparent_type,
	Light_Diffuse_type,
	Disney_type

};

void fixIndex(int& v, const int& n)
{
	v = v < 0 ? v + n : v > 0 ? v - 1 : -1;
}

void SkipSpace(char *&t)
{
	t += strspn(t, " \t");
}

/*struct NormalTexcoord
{
vec3 n[3];
vec3 t[3];
};*/

struct Light
{
	//pu vu em ne anh!
	vec3 p;
	vec3 u;
	vec3 v;
	vec3 n;
	vec3 e;
	vec3 area_pdf;
	//area  : area of this light source 
	//pdf   : the probability this light source will be chosen while sampling: = area / sum_area 
	////index : index of this light source in light_texs

	Light() {}
	Light(vec3 p_, vec3 u_, vec3 v_, vec3 e_, vec3 n_, vec3 area_pdf_) : p(p_), u(u_), v(v_), e(e_), n(n_), area_pdf(area_pdf_) {}
};

enum face
{
	Single_Line, Double_Line
};

void getdirection(string filepath, string& direction)
{
	size_t found = filepath.find_last_of('/\\');
	direction = filepath.substr(0, found + 1);
}

void getfilename(const string& filepath, string& filename)
{

	size_t found = filepath.find_last_of('/\\');
	filename = filepath.substr(found + 1);
}

void get_face_index(char*& t, const int& vs, const int& vts, const int& vns, vector<Face>& faces)
{
	string s = t;
	int length = s.find_last_of("0123456789");
	s = s.substr(0, length + 1);

	int sign = 1;
	int count = 0;
	vector<int> index;
	//vector<int> indices;
	int face_type = Single_Line;

	int num_data_per_vertex = 0;
	bool found_num_data_per_vertex = false;

	for (int i = 0; i <= length + 1; ++i)
	{
		if (s[i] == '-')
			sign = -1;
		else if (isdigit(s[i]))
		{
			count = 10 * count + s[i] - '0';
		}
		else if (s[i] == '/')
		{
			if (!found_num_data_per_vertex)
				++num_data_per_vertex;
			face_type = Single_Line;
			index.emplace_back(sign * count);
			sign = 1;
			count = 0;

			if (s[i + 1] == '/')
			{
				face_type = Double_Line;
				++i;
			}
		}
		else if (s[i] == ' ')
		{

			index.emplace_back(sign * count);
			sign = 1;
			count = 0;
			if (!found_num_data_per_vertex)
			{
				++num_data_per_vertex;
				found_num_data_per_vertex = true;
			}
		}
		else if (i == length + 1)
		{
			index.emplace_back(sign * count);
			sign = 1;
			break;
		}
	}


	int size = index.size();


	if (num_data_per_vertex == 3) // v/vt/vn case   12 18 24 30 
	{
		//cout << "line 3\n";
		for (int i = 0; i < size; i += 3)
		{
			fixIndex(index[i], vs);
			fixIndex(index[i + 1], vts);
			fixIndex(index[i + 2], vns);
		}

		int start_v = 0;
		int start_vt = 1;
		int start_vn = 2;

		int num_Triangle = size / 3 - 2;

		for (int i = 0; i < num_Triangle; ++i)
		{
			Face f(index[start_v], index[start_vt], index[start_vn], index[3 * i + 3], index[3 * i + 4], index[3 * i + 5], index[3 * i + 6], index[3 * i + 7], index[3 * i + 8]);
			faces.emplace_back(f);
		}
	}
	else if (num_data_per_vertex == 2)  //  v/vt or v//vn
	{
		//cout << "line 2\n";
		if (face_type == Single_Line) // v / vt
		{
			//cout << "single\n";
			for (int i = 0; i < size; i += 2)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vts);
			}

			int num_Triangle = size / 2 - 2;

			int start_v = 0;
			int start_vt = 1;

			for (int i = 0; i < num_Triangle; ++i)
			{
				//01 23 45    01 45 67
				Face f(index[start_v], index[start_vt], -1, index[2 * i + 2], index[2 * i + 3], -1, index[2 * i + 4], index[2 * i + 5], -1);
				faces.emplace_back(f);
			}
		}
		else if (face_type == Double_Line)// v // vn
		{
			//cout << "double\n";
			for (int i = 0; i < size; i += 2)
			{
				fixIndex(index[i], vs);
				fixIndex(index[i + 1], vns);
			}

			int num_Triangle = size / 2 - 2;

			int start_v = 0;
			int start_vn = 1;

			for (int i = 0; i < num_Triangle; ++i)
			{
				Face f(index[start_v], -1, index[start_vn], index[2 * i + 2], -1, index[2 * i + 3], index[2 * i + 4], -1, index[2 * i + 5]);
				faces.emplace_back(f);
			}
		}
	}
}

//convert image from size img_width, img_height to size width, height

//https://chao-ji.github.io/jekyll/update/2018/07/19/BilinearResize.html

void resize_image(unsigned char* original_image, int img_width, int img_height, vector<vec3>& resize_image, int width, int height)
{
	vector<vec3> image;

	float inv_255 = 1.0f / 255.0f;

	for (int i = 0; i < img_width * img_height * 3; i += 3)
	{
		float x = 255 * (original_image[i] * inv_255);
		float y = 255 * (original_image[i + 1] * inv_255);
		float z = 255 * (original_image[i + 2] * inv_255);

		image.emplace_back(vec3(x, y, z));
	}

	resize_image.resize(width * height);

	float x_ratio = width > 1 ? float(img_width) / (width) : 1;
	float y_ratio = height > 1 ? float(img_height) / (height) : 1;

	for (int i = 0; i < height; ++i)
	{
		for (int j = 0; j < width; ++j)
		{

			int xl = floor(x_ratio * j);
			int yl = floor(y_ratio * i);
			int xh = ceil(x_ratio * j);
			int yh = ceil(y_ratio * i);

			vec3 a = image[yl * img_width + xl];//yl, xl
			vec3 b = image[yl * img_width + xh];//yl, xh
			vec3 c = image[yh * img_width + xl];//yh, xl
			vec3 d = image[yh * img_width + xh];//yh, xh

			float x_weight = float(x_ratio * j - xl);
			float y_weight = float(y_ratio * i - yl);

			vec3 pixel = a * (1 - x_weight) * (1 - y_weight) +
				b * x_weight * (1 - y_weight) +
				c * y_weight * (1 - x_weight) +
				d * x_weight * y_weight;

			resize_image[i * width + j] = pixel;

		}
	}

	vector<vec3>().swap(image);
	delete[] original_image;
}

void copy_texture(vector<vec3>& source, vector<vec3>& destination)
{
	copy(source.begin(), source.end(), back_inserter(destination));

	vector<vec3>().swap(source);
}

struct Scene
{
	//---------------CPU data-------------

	int frame_count = 0;
	int leaf = 4;
	Camera camera;

	//Geometry
	vector<vec3> vertices;
	vector<vec3> normals;
	vector<vec2> texcoords;
	vector<Triangle> triangles;

	//transform model so the min bbox position will be at vec3(0, 0, 0)
	//this will make ray box intersection faster because many bbox will now have bbox value = 0 
	//vec3 transformation_vector;

	//Material
	vector<Material> mats;
	vector<Light> lights;

	//Texture
	//vector<Texture> textures;
	Texture tex;
	int require_width;
	int require_height;

	vector<unsigned char> albedo_textures_data;

	//unsigned char* albedo_textures_data;

	//vector<unsigned char> textures_data;

	unordered_map<string, int> texture_map;

	//BVH
	vector<FlatNode> flat_nodes;

	//Utility
	unordered_map<string, int> material_map;
	int material_count = 0;
	float texture_gammar = 2.6f;

	//----------------GPU data---------------
	GLuint vertices_tex, normals_tex, texcoords_tex, triangles_tex, mats_tex, lights_tex, bvh_tex;
	GLuint vertices_buffer, normals_buffer, texcoords_buffer, triangles_buffer, mats_buffer, lights_buffer, bvh_buffer;

	//Texture Handle
	GLuint albedo_textures;
	//texture
	//Texture textures;
	//GLuint albedo_texture, normal_texture, metalic_roughness_texture;


	//Shader
	Shader *path_trace_shader, *accumulate_shader, *output_shader, *path_trace_fade_shader;
	Quad* quad;

	//FBO
	GLuint path_trace_fbo, path_trace_texture;
	GLuint accumulate_fbo, accumulate_texture;
	GLuint path_trace_fade_fbo, path_trace_fade_texture;

	//bool half_res;

	Scene() {}
	Scene(string file_name, string shader_direction)
	{
		require_width = renderOptions.require_tex_width;
		require_height = renderOptions.require_tex_height;

		//vec3 v1();
		//vec3 v2(-0.02171, -0.0964, -0.9951);

		//vec3 sum = v1 + v2;

		//cout << sum.x << " " << sum.y << " " << sum.z << "\n";

		//vec3 look_from(-2.755610, 2.745992, 7.58545);
		//vec3 view_direction(0.0f, 0.001f, -1.0001f);
		//vec3 view_direction(0.0f, 0.0f, -1.0f);

		//camera = Camera(vec3(-2.70539, 2.49926, 2.15389), vec3(-2.7271, 2.40286, 1.15879), 40.0f);

		//camera = Camera(vec3(5.105184f, 0.731065f, -2.317890f), vec3(1.452592f, 1.013640f, -1.317287f), 70.0f);

		//cornell
		camera = Camera(vec3(-2.755610, 2.745992, 7.58545), vec3(-2.755610, 2.745992, 6.58545), 40.0f);
		
		//fire room
		//camera = Camera(vec3(5.105184f, 0.731065f, -2.317890f), vec3(1.452592f, 1.013640f, -1.317287f), 70.0f);

		//bath room
		//camera = Camera(vec3(0.0072405338287353516, 0.9424049544334412, -0.2275838851928711), vec3(-2.787562608718872, 0.9699121117591858, -2.6775901317596436), 55.0f);

		//Lamp
		//camera = Camera(vec3(7.755991f, 5.067980f, -6.643479f), vec3(0.238440f, 2.568779, 1.390990f), 36.0f);



		//vec3 look_from(0.0072405338287353516, 0.9424049544334412, -0.2275838851928711);
		//vec3 look_at(-2.787562608718872, 0.9699121117591858, -2.6775901317596436);
		//vec3 view_direction = normalize(look_at - look_from);
		//float fov = 55.0f;




		Read_Object(file_name);
		cout << "Complete Reading\n";

		build_bvh(leaf);
		cout << "Complte Building BVH\n";

		gpu_data(shader_direction);

		cout << vertices.size() << "\n";
		cout << normals.size() << "\n";
		cout << texcoords.size() << "\n";
		cout << triangles.size() << "\n";
		cout << flat_nodes.size() << "\n";

		delete_cpu_data();
		delete_gpu_data();
	}

	void ReadMtl(const string& filename, string& direction, unordered_map<string, int>& mtlMap, const string& EnviromentName = "")
	{
		ifstream f(filename);
		if (!f)
			cout << "Mtl file not exist\n";

		char line[256];
		int numberOfMaterial = 0;


		char tex_name[32];
		unordered_set<string> albedo_map;

		while (f.getline(line, 256))
		{
			if (line[0] == 'n')
				numberOfMaterial++;

			if (strncmp(line, "map_Kd", 6) == 0)
			{
				char* t = line;
				string real_name;
				sscanf_s(t += 7, "%s", tex_name);

				getfilename(string(tex_name), real_name);

				albedo_map.insert(real_name);
			}
		}

		tex.num_albedo = albedo_map.size();

		cout << "Texture size: " << tex.num_albedo << "\n";

		//tex.albedo.resize(require_width * require_height * 3 * tex.num_albedo);

		mats.resize(numberOfMaterial);

		f.clear();
		f.seekg(0, ios::beg);

		//int countTexture = -1;
		int countMaterial = -1;

		bool found_texture = false;

		char mtlname[64];
		char texturename[64];

		unordered_map<string, int> tex_map;

		int count_albedo = 0;
		int count_light = 0;

		while (f.getline(line, 256))
		{

			char* t = line;
			SkipSpace(t);

			if (strncmp(t, "newmtl", 6) == 0)
			{

				found_texture = false;

				sscanf_s(t += 7, "%s", mtlname);

				mtlMap[mtlname] = ++countMaterial;
			}
			else if (strncmp(t, "type", 4) == 0)
			{
				char type[64];
				sscanf_s(t += 5, "%s", type);
				if (strncmp(type, "Mirror", 6) == 0)
					mats[countMaterial].albedo.w = Mirror_type;
			}
			else if (t[0] == 'K')
			{
				
				vec4 emission(-1);
				if (t[1] == 'd')
					sscanf_s(t += 3, "%f %f %f", &mats[countMaterial].albedo.x, &mats[countMaterial].albedo.y, &mats[countMaterial].albedo.z);
				if (t[1] == 'e')
					sscanf_s(t += 3, "%f %f %f", &emission.x, &emission.y, &emission.z);

				if (emission.x > 0 || emission.y > 0 || emission.z > 0)
					emission.w = count_light++;

				mats[countMaterial].emission = emission;
			}
			else if (strncmp(t, "map_Kd", 6) == 0)
			{
				string real_name;
				sscanf_s(t += 7, "%s", texturename);

				getfilename(string(texturename), real_name);

				if (!texture_map.count(real_name))
				{
					//cout << real_name << "\n";

					texture_map[real_name] = count_albedo;
					mats[countMaterial].tex_ind.x = count_albedo;

					//cout << "albedo index: " << count_albedo << "\n";

					string name = direction + real_name;

					//cout << name << "\n";

					int w, h;

					unsigned char* texture = stbi_load(name.c_str(), &w, &h, 0, 3);


					//int tex_byte = require_width * require_height * 3;

					tex.albedo_size = vec2(require_width, require_height);

					//unsigned char* resize_texture = new unsigned char[tex_byte];

					//cout << w << " " << h << " " << require_width << " " << require_height << "\n";

					vector<vec3> resize_texture;// (require_width * require_height);
					if (w != require_width || h != require_height)
					{
						//resize_image(texture, w, h, resize_texture, require_width, require_height);
						//(&texture[0], w, h, 0, resize_texture, require_width, require_height, 0, 3);

						//copy(resize_texture.begin(), resize_texture.end(), &tex.albedo[count_albedo * tex_byte]);
						//cout << "Resize Success";

						resize_image(texture, w, h, resize_texture, require_width, require_height);

						for (int i = 0; i < require_width * require_height; ++i)
						{
							vec3 v = resize_texture[i];
							//v = vec3(powf(v.x, 1 / 2.2f), powf(v.y, 1 / 2.2f), powf(v.z, 1 / 2.2f));
							//v = vec3(powf(v.x, 2.2f), powf(v.y, 2.2f), powf(v.z, 2.2f));
							tex.albedo.emplace_back(v);
							//cout << v.x << " " << v.y << " " << v.z << "\n";
						}
						//copy_texture(resize_texture, tex.albedo);
						//cout << "Pass Copy\n";
					}
					else
					{
						float inv_255 = 1.0f / 255.0f;
						for (int i = 0; i < 3 * require_width * require_height; i += 3)
						{
							float x = 255 * (texture[i] * inv_255);
							float y = 255 * (texture[i + 1] * inv_255);
							float z = 255 * (texture[i + 2] * inv_255);

							//cout << x << " " << y << " " << z << "\n";

							//vec3 v = vec3(powf(x, 2.2f), powf(y, 2.2f), powf(z, 2.2f));
							//tex.albedo.emplace_back(v);

							tex.albedo.emplace_back(vec3(x, y, z));
						}
					}

					//else
					//	copy()
					//getchar();
					//memcpy(&(tex.albedo[require_width * require_height * 3 * count_albedo]), &resize_texture[0], 3 * require_width * require_height);
					//memcpy(resize_texture, resize_texture + tex_byte, &tex.albedo[count_albedo * tex_byte]);
					//delete[] texture;

					++count_albedo;

					//cout << "pass\n";
					//else
					//	copy(tex.albedo.begin(), tex.albedo.end(), &)
				}
			}
		}

		//copy everything to albedo_textures_data

		if (tex.num_albedo > 0)
		{
			//int s = 0;
			//albedo_textures_data = new unsigned char[3 * require_width * require_height * tex.num_albedo]();

			for (auto& v : tex.albedo)
			{
				float x = v.x;
				float y = v.y;
				float z = v.z;

				//if (x != 0 || y != 0 || z != 0)
				//cout << x << " " << y << " " << z << "\n";

				//albedo_textures_data[s] = (unsigned char)(x * 255);
				//albedo_textures_data[s + 1] = (unsigned char)(y * 255);
				//albedo_textures_data[s + 2] = (unsigned char)(z * 255);
				//s += 3;
				albedo_textures_data.emplace_back(x);
				albedo_textures_data.emplace_back(y);
				albedo_textures_data.emplace_back(z);
			}
		}

		vector<vec3>().swap(tex.albedo);

		//getchar();
		//Resize Texture


		//int require_width = renderOptions.require_tex_width;
		//int require_height = renderOptions.require_tex_height;
		//int texture_bytes = require_width * require_height * 4;


		//textures_data.resize(texture_bytes * textures.size());

		/*for (int i = 0; i < textures.size(); ++i)
		{
		int tex_width = textures[i].w;
		int tex_height = textures[i].h;

		if (tex_width != renderOptions.require_tex_width || tex_height != renderOptions.require_tex_height)
		{
		unsigned char* resize_texture = new unsigned char[texture_bytes];
		stbir_resize_uint8(&textures[i].tex_data[0], tex_width, tex_height, 0, resize_texture, require_width, require_height, 0, 4);
		copy(resize_texture, resize_texture + texture_bytes, &textures_data);
		delete[] resize_texture;
		}
		else
		copy(textures[i].tex_data.begin(), textures[i].tex_data.end(), &textures_data[i * texture_bytes]);
		}*/
	}

	void Read_Object(string file_name, string enviroment_path = "")
	{
		//cout << file_name << "\n";
		ifstream f(file_name);
		if (!f)
			cout << "obj file not exist";

		string direction;
		getdirection(file_name, direction);

		//cout << direction << "\n";

		char line[1024], mtl_name[256];
		char mtlib[64];

		int num_v = 0, num_vt = 0, num_vn = 0;

		int mtl_ind;

		unordered_map<string, int> mtlMap;

		bool read_mtl = false;

		int light_index = 0;

		vec3 vertex_min(inf, inf, inf);
		vec3 vertex_max(-inf, -inf, -inf);

		while (f.getline(line, 1024))
		{
			char* t = line;
			int prev_space = strspn(t, " \t");
			t += prev_space;

			if (strncmp(t, "v", 1) == 0)
			{
				t += 1;
				float x, y, z;
				if (strncmp(t, " ", 1) == 0)
				{
					int post_space = strspn(t, " \t");
					t += post_space;
					sscanf_s(t, "%f %f %f", &x, &y, &z);

					vec3 vertex(x, y, z);

					vertex_min = min_vec(vertex_min, vertex);
					vertex_max = max_vec(vertex_max, vertex);

					vertices.emplace_back(vertex);
					++num_v;
				}
				else if (strncmp(t, "t", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					sscanf_s(t, "%f %f", &x, &y);
					texcoords.emplace_back(vec2(x, 1.0f - y));
					++num_vt;
				}
				else if (strncmp(t, "n", 1) == 0)
				{
					t += 1;
					int post_space = strspn(t, " \t");
					t += post_space;

					sscanf_s(t, "%f %f %f", &x, &y, &z);
					normals.emplace_back(vec3(x, y, z));
					++num_vn;
				}
			}
			else if (strncmp(t, "f", 1) == 0)
			{
				//t += 1;
				int post_space = strspn(t + 1, " \t");

				t += post_space + 1;

				int face_type;
				vector<Face> faces;
				get_face_index(t, num_v, num_vt, num_vn, faces);

				for (int i = 0; i < faces.size(); ++i)
				{
					int v0 = faces[i].v[0];
					int v1 = faces[i].v[1];
					int v2 = faces[i].v[2];

					int vt0 = faces[i].vt[0];
					int vt1 = faces[i].vt[1];
					int vt2 = faces[i].vt[2];

					int vn0 = faces[i].vn[0];
					int vn1 = faces[i].vn[1];
					int vn2 = faces[i].vn[2];

					//if not exist, compute

					Triangle tr;
					if (vn0 == -1)
					{
						vec3 p0 = vertices[v0];
						vec3 p1 = vertices[v1];
						vec3 p2 = vertices[v2];

						vec3 n = cross((p1 - p0), (p2 - p0));

						//tr = Triangle(ivec4(v0, v1, v2, mtl_ind), ivec4(vt0, vt1, vt2, 0), ivec4(n.x, n.y, n.z, 0));
						tr = Triangle(ivec4(v0, v1, v2, mtl_ind), ivec4(n.x, n.y, n.z, 0), ivec4(vt0, vt1, vt2, 0));
					}
					else
						tr = Triangle(ivec4(v0, v1, v2, mtl_ind), ivec4(vn0, vn1, vn2, 1), ivec4(vt0, vt1, vt2, 0));
					if (mats[mtl_ind].emission.w != -1)//light
					{
						vec3 p0 = vertices[v0];
						vec3 p1 = vertices[v1];
						vec3 p2 = vertices[v2];

						vec3 u = p1 - p0;
						vec3 v = p2 - p0;
						vec3 e = vec3(mats[mtl_ind].emission.x, mats[mtl_ind].emission.y, mats[mtl_ind].emission.z);
						vec3 n = cross(u, v);

						//cout << e.x << " " << e.y << " " << e.z << "\n";

						float area = length(n);

						n = normalize(n);

						//tr.vt.w = mats[mtl_ind].emission.w;//light_index++;

						Light light(p0, u, v, e, n, vec3(area, 0, 0));

						lights.emplace_back(light);
					}

					triangles.emplace_back(tr);
				}
			}
			else if (strncmp(t, "usemtl", 6) == 0)
			{
				sscanf_s(t += 7, "%s", mtl_name);
				string realname = (string)mtl_name;

				mtl_ind = mtlMap[realname];
			}
			else if (t[0] == 'm' && !read_mtl)
			{
				sscanf_s(t += 7, "%s", mtlib);

				string realname = direction + (string)mtlib;


				ReadMtl(realname, direction, mtlMap, enviroment_path);
				read_mtl = true;
			}
		}

		//compute light pdf;

		float sum_area = 0.0f;

		for (auto& v : lights)
			sum_area += v.area_pdf.x;
		if (sum_area > 0.0f)
		{
			float inv_sum_area = 1.0f / sum_area;
			for (auto& v : lights)
				v.area_pdf.y = v.area_pdf.x * inv_sum_area;
		}

		vec3 center = (vertex_min + vertex_max) * 0.5f;

		vec3 transformation_vector = -vertex_min;//center - vertex_min;
		//now everything have to be update
		for (auto& v : vertices)
			v += transformation_vector;
		for (auto& v : lights)
			v.p += transformation_vector;

		camera.position += transformation_vector;
		cout << "Complete Transform\n";
	}

	//BFS build is equal to DFS build in term of speed
	void build_bvh(int leaf)
	{
		cout << "Start Buiding BVH\n";
		
		SBVH sbvh(triangles, vertices);

		cout << "Finish Building BVH\n";

		//cout << "Start Optimize!\n";

		//optimize(sbvh);

		//cout << "End Optimize\n";

		flat_nodes = sbvh.flat_nodes;


		/*vector<Triangle> new_trs;

		int s = sbvh.triangle_indices.size();

		new_trs.resize(s);

		for (int i = 0; i < s; ++i)
			new_trs[i] = triangles[sbvh.triangle_indices[i]];

		triangles = new_trs;

		vector<Triangle>().swap(new_trs);*/
		sbvh.clean_up();
	}

	void delete_cpu_data()
	{
		vector<vec3>().swap(vertices);
		vector<vec3>().swap(normals);
		vector<vec2>().swap(texcoords);

		vector<Triangle>().swap(triangles);

		vector<Material>().swap(mats);

		vector<Light>().swap(lights);

		vector<FlatNode>().swap(flat_nodes);

		unordered_map<string, int>().swap(material_map);
	}

	void delete_gpu_data()
	{
		glDeleteBuffers(1, &vertices_buffer);
		glDeleteBuffers(1, &normals_buffer);
		glDeleteBuffers(1, &texcoords_buffer);
		glDeleteBuffers(1, &triangles_buffer);
		glDeleteBuffers(1, &mats_buffer);
		glDeleteBuffers(1, &lights_buffer);
		glDeleteBuffers(1, &bvh_buffer);
	}

	void delete_tex_data()
	{
		glDeleteTextures(1, &vertices_tex);
		glDeleteTextures(1, &normals_tex);
		glDeleteTextures(1, &texcoords_tex);
		glDeleteTextures(1, &triangles_tex);
		glDeleteTextures(1, &mats_tex);
		glDeleteTextures(1, &lights_tex);
		glDeleteTextures(1, &bvh_tex);
	}

	void gpu_data(string shader_direction)
	{

		GLint maxAttach = 0;
		glGetIntegerv(GL_MAX_COLOR_ATTACHMENTS, &maxAttach);

		GLint maxDrawBuf = 0;
		glGetIntegerv(GL_MAX_DRAW_BUFFERS, &maxDrawBuf);

		cout << "max Attach:" << maxAttach << "\n";
		cout << "max Draw:" << maxDrawBuf << "\n";

		//GLuint vertices_tex, normal_texcoords_tex, triangles_tex, mats_tex, lights_tex;
		//GLuint vertices_buffer, normal_texcoords_buffer, triangles_buffer, mats_buffer, lights_buffer;

		glGenBuffers(1, &vertices_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, vertices_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(vec3) * vertices.size(), &vertices[0], GL_STATIC_DRAW);
		glGenTextures(1, &vertices_tex);
		glBindTexture(GL_TEXTURE_BUFFER, vertices_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, vertices_buffer);

		glGenBuffers(1, &normals_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, normals_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(vec3) * normals.size(), &normals[0], GL_STATIC_DRAW);
		glGenTextures(1, &normals_tex);
		glBindTexture(GL_TEXTURE_BUFFER, normals_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, normals_buffer);

		glGenBuffers(1, &texcoords_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, texcoords_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(vec2) * texcoords.size(), &texcoords[0], GL_STATIC_DRAW);
		glGenTextures(1, &texcoords_tex);
		glBindTexture(GL_TEXTURE_BUFFER, texcoords_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RG32F, texcoords_buffer);

		glGenBuffers(1, &triangles_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, triangles_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(Triangle) * triangles.size(), &triangles[0], GL_STATIC_DRAW);
		glGenTextures(1, &triangles_tex);
		glBindTexture(GL_TEXTURE_BUFFER, triangles_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32I, triangles_buffer);

		glGenBuffers(1, &mats_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, mats_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(Material) * mats.size(), &mats[0], GL_STATIC_DRAW);
		glGenTextures(1, &mats_tex);
		glBindTexture(GL_TEXTURE_BUFFER, mats_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, mats_buffer);

		glGenBuffers(1, &lights_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, lights_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(Light) * lights.size(), &lights[0], GL_STATIC_DRAW);
		glGenTextures(1, &lights_tex);
		glBindTexture(GL_TEXTURE_BUFFER, lights_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGB32F, lights_buffer);

		glGenBuffers(1, &bvh_buffer);
		glBindBuffer(GL_TEXTURE_BUFFER, bvh_buffer);
		glBufferData(GL_TEXTURE_BUFFER, sizeof(FlatNode) * flat_nodes.size(), &flat_nodes[0], GL_STATIC_DRAW);
		glGenTextures(1, &bvh_tex);
		glBindTexture(GL_TEXTURE_BUFFER, bvh_tex);
		glTexBuffer(GL_TEXTURE_BUFFER, GL_RGBA32F, bvh_buffer);

		//texture data
		if (tex.num_albedo > 0)
		{
			cout << "init albedo texture\n";
			glGenTextures(1, &albedo_textures);
			glActiveTexture(GL_TEXTURE0);
			glBindTexture(GL_TEXTURE_2D_ARRAY, albedo_textures);
			//tex.albedo_size.x, .y
			glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGB8, require_width, require_height, tex.num_albedo, 0, GL_RGB, GL_UNSIGNED_BYTE, &albedo_textures_data[0]);//&albedo_textures_data[0]);//&tex.albedo[0]);
																																								//glTexImage3D(GL_TEXTURE_2D_ARRAY, 0, GL_RGB, tex.albedo_size.x, tex.albedo_size.y, tex.num_albedo, 0, GL_RGB, GL_FLOAT, &tex.albedo[0]);//&tex.albedo[0]);
			glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
			glTexParameteri(GL_TEXTURE_2D_ARRAY, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
			glBindTexture(GL_TEXTURE_2D_ARRAY, 0);
			cout << "finish init albedo \n";
		}

		path_trace_shader = new Shader(shader_direction + "path_trace.vs", shader_direction + "path_trace.fs");
		//path_trace_shader = new Shader(shader_direction + "path_trace.vs", shader_direction + "path_trace_optimize_3_bbox.fs");
		accumulate_shader = new Shader(shader_direction + "accumulate.vs", shader_direction + "accumulate.fs");
		output_shader = new Shader(shader_direction + "output.vs", shader_direction + "output.fs");
		//path_trace_fade_shader = new Shader(shader_direction + "path_trace_fade.vs", shader_direction + "path_trace_fade.fs");

		quad = new Quad();

		//path trace fbo
		glGenFramebuffers(1, &path_trace_fbo);
		glBindFramebuffer(GL_FRAMEBUFFER, path_trace_fbo);

		//Create Texture for FBO
		glGenTextures(1, &path_trace_texture);
		glBindTexture(GL_TEXTURE_2D, path_trace_texture);
		cout << "Tex " << path_trace_texture << "\n";
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, screen_width, screen_height, 0, GL_RGB, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, path_trace_texture, 0);

		//path trace fade fbo
		glGenFramebuffers(1, &path_trace_fade_fbo);
		glBindFramebuffer(GL_FRAMEBUFFER, path_trace_fade_fbo);

		//create half texture
		glGenTextures(1, &path_trace_fade_texture);
		glBindTexture(GL_TEXTURE_2D, path_trace_fade_texture);
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, screen_width / 2, screen_height / 2, 0, GL_RGB, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, path_trace_fade_texture, 0);

		//accumulate fbo
		glGenFramebuffers(1, &accumulate_fbo);
		glBindFramebuffer(GL_FRAMEBUFFER, accumulate_fbo);

		//Create Texture for accumulate FBO
		glGenTextures(1, &accumulate_texture);
		glBindTexture(GL_TEXTURE_2D, accumulate_texture);

		cout << "Tex " << accumulate_texture << "\n";

		//int x = 0;
		glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB32F, screen_width, screen_height, 0, GL_RGB, GL_FLOAT, 0);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
		glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
		glBindTexture(GL_TEXTURE_2D, 0);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, accumulate_texture, 0);

		path_trace_shader->use();
		//path_trace_shader->setInt("accumulate_tex", 0);
		path_trace_shader->setInt("vertices_tex", 1);
		path_trace_shader->setInt("normals_tex", 2);
		path_trace_shader->setInt("texcoords_tex", 3);
		path_trace_shader->setInt("triangles_tex", 4);
		path_trace_shader->setInt("mats_tex", 5);
		path_trace_shader->setInt("lights_tex", 6);
		path_trace_shader->setInt("bvh", 7);
		path_trace_shader->setInt("albedo_textures", 8);

		path_trace_shader->setVec3("camera.up", camera.up);
		path_trace_shader->setVec3("camera.right", camera.right);
		path_trace_shader->setVec3("camera.forward", camera.forward);
		path_trace_shader->setVec3("camera.position", camera.position);
		path_trace_shader->setFloat("camera.fov", camera.fov);
		path_trace_shader->setFloat("camera.focalDist", camera.focalDist);
		path_trace_shader->setFloat("camera.aperture", camera.aperture);

		path_trace_shader->setVec2("screenResolution", vec2(screen_width, screen_height));
		path_trace_shader->setInt("numLights", lights.size());
		path_trace_shader->setInt("maxDepth", 3);
		path_trace_shader->reset();

	}

	void Render()
	{
		if (camera.isMoving)
		{
			glBindFramebuffer(GL_FRAMEBUFFER, path_trace_fbo);
			glClear(GL_COLOR_BUFFER_BIT);

			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			glBindFramebuffer(GL_FRAMEBUFFER, accumulate_fbo);
			glClear(GL_COLOR_BUFFER_BIT);

			glBindFramebuffer(GL_FRAMEBUFFER, 0);

			frame_count = 0;
			//bo vo la failed!
			//glActiveTexture(GL_TEXTURE0);
			//glBindTexture(GL_TEXTURE_2D, accumulate_texture);
			//glClear(GL_COLOR_BUFFER_BIT);
		}

		//accumulate chi co vai tro luu tam thoi cho 1 frame
		//de no o day thi hinh se bi lo ro

		//un comment for single sample or 1 sps rendering!
		/*
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, accumulate_texture);
		glClear(GL_COLOR_BUFFER_BIT);
		*/
		glActiveTexture(GL_TEXTURE1);
		glBindTexture(GL_TEXTURE_BUFFER, vertices_tex);
		glActiveTexture(GL_TEXTURE2);
		glBindTexture(GL_TEXTURE_BUFFER, normals_tex);
		glActiveTexture(GL_TEXTURE3);
		glBindTexture(GL_TEXTURE_BUFFER, texcoords_tex);
		glActiveTexture(GL_TEXTURE4);
		glBindTexture(GL_TEXTURE_BUFFER, triangles_tex);
		glActiveTexture(GL_TEXTURE5);
		glBindTexture(GL_TEXTURE_BUFFER, mats_tex);
		glActiveTexture(GL_TEXTURE6);
		glBindTexture(GL_TEXTURE_BUFFER, lights_tex);
		glActiveTexture(GL_TEXTURE7);
		glBindTexture(GL_TEXTURE_BUFFER, bvh_tex);
		glActiveTexture(GL_TEXTURE8);
		glBindTexture(GL_TEXTURE_2D_ARRAY, albedo_textures);
		float r1 = 0, r2 = 0, r3 = 0;
		//r1 = ((float)rand() / (RAND_MAX)), r2 = ((float)rand() / (RAND_MAX)), r3 = ((float)rand() / (RAND_MAX));
		//r1 = randf(), r2 = randf(), r3 = randf();

		r1 = randf2(), r2 = randf2();// , r3 = randf2();
		path_trace_shader->use();
		path_trace_shader->setVec2("randomVector", vec2(r1, r2));// , r3));

		glBindFramebuffer(GL_FRAMEBUFFER, path_trace_fbo);
		quad->Draw(path_trace_shader);

		//glBindFramebuffer(GL_FRAMEBUFFER, 0);

		glBindFramebuffer(GL_FRAMEBUFFER, accumulate_fbo);
		quad->Draw(accumulate_shader);

		glBindFramebuffer(GL_FRAMEBUFFER, 0);


		++frame_count;
		float inv_frame_count = (float)1.0f / frame_count;

		output_shader->use();
		output_shader->setFloat("invSampleCounter", inv_frame_count);
		glActiveTexture(GL_TEXTURE0);
		glBindTexture(GL_TEXTURE_2D, path_trace_texture);
		quad->Draw(output_shader);
	}

	void update(float second)
	{
		path_trace_shader->use();

		path_trace_shader->setVec3("camera.position", camera.position);
		path_trace_shader->setVec3("camera.right", camera.right);
		path_trace_shader->setVec3("camera.up", camera.up);
		path_trace_shader->setVec3("camera.forward", camera.forward);
		path_trace_shader->setFloat("camera.fov", camera.fov);
		path_trace_shader->setFloat("camera.focalDist", camera.focalDist);
		path_trace_shader->setFloat("camera.aperture", camera.aperture);

		path_trace_shader->reset();
	}


};

#endif // !_SCENE_SANG_H_

