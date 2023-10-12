#version 460

in vec2 tex;

out vec3 color;
uniform sampler2D accumulate_tex;
uniform samplerBuffer vertices_tex;
uniform samplerBuffer normals_tex;
uniform samplerBuffer texcoords_tex;
uniform isamplerBuffer triangles_tex; //isamplerBuffer: reduce float to int conversion, result in 2.8% speed up!
uniform samplerBuffer mats_tex; 
uniform samplerBuffer lights_tex;
uniform samplerBuffer bvh;
uniform sampler2DArray albedo_textures;

#define pi 3.1415926
#define pi2 6.2831853
#define ipi 1.0f / pi

struct Ray { vec3 o; vec3 d; };
struct Camera { vec3 up; vec3 right; vec3 forward; vec3 position; float fov; float focalDist; float aperture; };


struct Light { vec3 position; vec3 emission; vec3 u; vec3 v; vec3 area; };


vec2 seed;

uniform Camera camera;
uniform vec2 screenResolution;
uniform vec2 randomVector;

uniform int numLights;

#define inf 1e9
#define eps 1e-4f

float rand()
{
	seed -= vec2(randomVector.x * randomVector.y);
	return fract(sin(dot(seed, vec2(12.9898, 78.233))) * 43758.5453);
}

void onb(in vec3 n, inout vec3 u, inout vec3 v)
{	
	if (n.z < -0.9999999f) // Handle the singularity
	{
		u = vec3(0.0f, -1.0f, 0.0f);
		v = vec3(-1.0f, 0.0f, 0.0f);
		return;
	}
	else
	{
		float a = 1.0f / (1.0f + n.z);
		float b = -n.x * n.y * a;
		u = vec3(1.0f + b, b, -n.x);
		v = vec3(b, 1.0f + b, -n.y);
		return;
	}
}

void min_max(inout vec3 bmin, inout vec3 bmax)
{
	if(bmin.x > bmax.x)
	{
		float t = bmax.x;
		bmax.x = bmin.x;
		bmin.x = t;
	}
	if(bmin.y > bmax.y)
	{
		float t = bmax.y;
		bmax.y = bmin.y;
		bmin.y = t;
	}
	if(bmin.z > bmax.z)
	{
		float t = bmax.z;
		bmax.z = bmin.z;
		bmin.z = t;
	}
}

float hit_bbox(Ray r, vec3 bmin, vec3 bmax, vec3 invdir, out float tl)
{
	
	/*bmax[0] = (bmax[0] - r.o[0]) * invdir[0];
	bmax[1] = (bmax[1] - r.o[1]) * invdir[1];
	bmax[2] = (bmax[2] - r.o[2]) * invdir[2];

	bmin[0] = (bmin[0] - r.o[0]) * invdir[0];
	bmin[1] = (bmin[1] - r.o[1]) * invdir[1];
	bmin[2] = (bmin[2] - r.o[2]) * invdir[2];*/


	bmin = (bmin - r.o) * invdir;
	bmax = (bmax - r.o) * invdir;
		
	//vec3 tmin = min(bmin, bmax), tmax = max(bmin, bmax);
	
	vec3 tmax = max(bmax, bmin);
	bmin = min(bmax, bmin);

	float th = min(tmax.x, min(tmax.y, tmax.z));
	
	tl = max(bmin.x, max(bmin.y, bmin.z));

	return th;
}

/*void tmin_tmax(in vec3 f, in vec3 n, out vec3 tmax, out vec3 tmin)
{
	tmax = f;
	tmin = n;

	if(f.x < n.x)
	{
		tmax.x = n.x;
		tmin.x = f.x;
	}
	if(f.y < n.y)
	{
		tmax.y = n.y;
		tmin.y = f.y;
	}
	if(f.z < n.z)
	{
		tmax.z = n.z;
		tmin.z = f.z;
	}
}

float hit_bbox(in vec3 o, in vec3 bmin, in vec3 bmax, in vec3 invdir, out float tl)
{
	//vec3 invdir = 1.0f / r.d;

	//vec3 f = (bmax - r.o) * invdir, n = (bmin - r.o) * invdir;

	bmin = (bmin - o) * invdir, bmax = (bmax - o) * invdir;


	vec3 tmax = max(bmax, bmin), tmin = min(bmax, bmin);

	
	tl = max(tmin.x, max(tmin.y, tmin.z));

	return min(tmax.x, min(tmax.y, tmax.z));
}
*/



/*float hit_bbox(Ray r, vec3 bmin, vec3 bmax, vec3 invdir, out float tl, float max_t)
{
	//vec3 invdir = 1.0f / r.d;

	vec3 f = (bmax - r.o) * invdir, n = (bmin - r.o) * invdir;

	vec3 tmax = max(f, n), tmin = min(f, n);

	//vec3 tmax, tmin;

	//tmin_tmax(f, n, tmax, tmin);

	//float th = min(tmax.x, min(tmax.y, tmax.z));
	tl = max(tmin.x, max(tmin.y, tmin.z));

	return tl < max_t ? min(tmax.x, min(tmax.y, tmax.z)) : -1.0f;

	//if(tl > max_t)
	//	return -1;
	//return th >= 0.0f && (tl <= th);

	//return (1.00000024f * t1 >= t0) ? (t0 >= eps ? t0 : t1) : -1.0;

	//return min(tmax.x, min(tmax.y, tmax.z));
	//return th;
	//return (th >= tl) ? (tl >= 0.0f ? tl : th ) : -1.0;
}
*/
/*struct Material
{
	vec3 albedo;
	vec3 emission;
};

struct LightSampleRec 
{ 
	vec3 position; 
	vec3 normal; 
	vec3 emission; 
	float pdf; 
};
struct BSDFSampleRec
{
	vec3 direction;
	float pdf;
	bool is_specular;
};


void SampleQuad(in Light light, inout LightSampleRec light_rec)
{
	float r1 = rand();
	float r2 = rand();

	light_rec.position = light.position + light.u * r1 + light.v * r2;
	light_rec.normal = normalize(cross(light.u, light.v));
	light_rec.emission = light.emission * numLights;
}
*/


float power_heuristic(float a, float b)
{
	float t = a * a;
	return t / (b * b + t);
}

/*struct Material
{
	vec3 albedo;
	vec3 emission;
}*/

/*struct Material
{
	vec4 albedo;  //xyz, w = alpha value, 0 = complete transparent, 1 = complete opaque
	vec4 emission;//xyz, w = -1 : not emission, w = 0, 1, .... : emission and light_index
	vec4 specular;//xyz, w = 0 : not specular matrial, w = 1 : specular material
	vec4 tex_ind;//albedo ind, normal_ind, specular_ind, metallic_roughness_ind
};*/	

struct HitRecord
{
	//Material mat;
	vec4 albedo;  //xyz, w = alpha value, 0 = complete transparent, 1 = complete opaque
	vec4 emission;//xyz, w = -1 : not emission, w = 0, 1, .... : emission and light_index
	vec4 specular;//xyz, w = 0 : not specular matrial, w = 1 : specular material
	vec4 tex_ind;//albedo ind, normal_ind, specular_ind, metallic_roughness_ind

	vec3 n;
	vec2 texcoord;
	float u;
	
	float v;
	float t;
	int triangle_ind;
	//float to int conversion is really expensive
	//so we decleare mtl_ind as float to avoid many conversion during ray triangle intersection
	//we only need to convert once, after bvh traversal complete
	int mtl_ind;
	//int mtl_ind;
};


vec3 cosine_hemisphere_sampling()
{
	float u1 = rand();
	float u2 = rand();
	
	float r = sqrt(u1);
	float phi = pi2 * u2;

	//x = r * cos(phi);
	//y = r * sin(phi)
	//z = sqrt(1 - x * x - y * y) = sqrt(1 - r ^ 2) = sqrt(1 - u1)

	return vec3(r * cos(phi), r * sin(phi), sqrt(1.0f - u1));
}

//---------Diffuse Lambert-------------

vec3 diffuse_sample(in Ray r, vec3 n)//inout HitRecord rec)
{
	//vec3 n = rec.n;

	//vec3 up = abs(n.z) < 0.999 ? vec3(0, 0, 1) : vec3(1, 0, 0);
	//vec3 u = normalize(cross(up, n));
	//vec3 v = cross(n, u);

	vec3 u, v;
	onb(n, u, v);

	vec3 dir = cosine_hemisphere_sampling();
	dir = u * dir.x + v * dir.y + n * dir.z;

	return dir;
}

float diffuse_pdf(Ray r, inout HitRecord rec, in vec3 sample_direction)
{
	float cos = dot(sample_direction, rec.n);
	return cos * ipi;

	//return cos > eps ? cos * ipi : 0.0f;
}

vec3 diffuse_bsdf(Ray r, inout HitRecord rec, in vec3 sample_direction)
{
	//float cos_o = dot(sample_direction, rec.n);

	//int mtl_ind = rec.mtl_ind;
	//vec3 albedo = rec.mat.albedo.xyz;//texelFetch(mats_tex, 3 * mtl_ind).xyz;

	vec3 albedo = rec.albedo.xyz;

	return albedo;
	//return cos_o > eps ? cos_o * albedo  : vec3(0.0f);
} 

vec2 interpolate(vec2 a, vec2 b, vec2 c, float u, float v)
{
	return a * (1.0f - u - v) + b * u + c * v;
}

vec3 interpolate(vec3 a, vec3 b, vec3 c, float u, float v)
{
	return a * (1.0f - u - v) + b * u + c * v;
}

bool hit_triangle(Ray r, int i, inout HitRecord rec, inout float t)
{
	ivec4 triIndex = texelFetch(triangles_tex, 3 * i);
	
	vec3 v0 = texelFetch(vertices_tex, triIndex.x).xyz; // int(triIndex.x)
	vec3 v1 = texelFetch(vertices_tex, triIndex.y).xyz;
	vec3 v2 = texelFetch(vertices_tex, triIndex.z).xyz;

	//replace
	//e0 = v1 
	//e1 = v2;

	//vec3 e0 = v1 - v0;
	//vec3 e1 = v2 - v0;
	
	v1 -= v0;
	v2 -= v0;

	vec3 pv = cross(r.d, v2);
	

	vec3 tv = r.o - v0;
	vec3 qv = cross(tv, v1);

	vec4 uvt;
	uvt.x = dot(tv, pv);//u
	//if(uvt.x < 0.0f)
	//	return false;

	uvt.y = dot(r.d, qv);//v
	//if(uvt.y < 0.0f)
	//	return false;

	uvt.z = dot(v2, qv);//t
	//float det = dot(e0, pv);
	float inv_det = 1.0f / dot(v1, pv);//det;

	uvt.xyz = uvt.xyz * inv_det;
	uvt.w = 1.0 - uvt.x - uvt.y;//u + v >= det

	
	if (all(greaterThanEqual(uvt, vec4(0.0f))) && uvt.z < t)
	{
		t = uvt.z;
		//rec.t = t;
		rec.u = uvt.x;
		rec.v = uvt.y;
		rec.triangle_ind = i;
		rec.mtl_ind = (triIndex.w);
		return true;
	}	
	return false;
}

bool hit_triangle_no_rec(Ray r, int i, float max_t)
{
	ivec4 triIndex = texelFetch(triangles_tex, 3 * i);
	
	vec3 v0 = texelFetch(vertices_tex, triIndex.x).xyz; // int(triIndex.x)
	vec3 v1 = texelFetch(vertices_tex, triIndex.y).xyz;
	vec3 v2 = texelFetch(vertices_tex, triIndex.z).xyz;

	//vec3 e0 = v1 - v0;
	//vec3 e1 = v2 - v0;
	v1 -= v0;
	v2 -= v0;

	vec3 pv = cross(r.d, v2);
	//float det = dot(e0, pv);

	vec3 tv = r.o - v0;
	vec3 qv = cross(tv, v1);

	vec4 uvt;
	uvt.x = dot(tv, pv);//u
	//if(uvt.x < 0.0f)
	//	return false;

	uvt.y = dot(r.d, qv);//v
	//if(uvt.y < 0.0f)
	//	return false;

	uvt.z = dot(v2, qv);//t
	float inv_det = 1.0f / dot(v1, pv);//1.0f / det;

	uvt.xyz = uvt.xyz * inv_det;
	uvt.w = 1.0 - uvt.x - uvt.y;//u + v >= det

		
	return (uvt.z < max_t && all(greaterThanEqual(uvt, vec4(0.0f))));
}

void compute_normal_texcoord_and_material(inout HitRecord rec)
{
	//triangle ind alow access to vn and vt
	//mtl_ind allow acess to properties of material this triangle hold

	int triangle_ind = rec.triangle_ind;

	int mtl_ind = (rec.mtl_ind);

	//vec4 tex_ind = texelFetch(mats_tex, 4 * mtl_ind + 3);

	

	//no texture
	//rec.texcoord.x = -1;
	
	//use texture
	/*if(tex_ind.x != -1)
	{
		vec2 vt0 = texelFetch(texcoords_tex, int(vt.x)).xy;
		vec2 vt1 = texelFetch(texcoords_tex, int(vt.y)).xy;
		vec2 vt2 = texelFetch(texcoords_tex, int(vt.z)).xy;

		rec.texcoord = interpolate(vt0, vt1, vt2, rec.u, rec.v);
	}*/
	
	ivec4 vn = texelFetch(triangles_tex, 3 * triangle_ind + 1);

	
	//0 mean no need to interpolate, just use directly
	if(vn.w == 0)
		rec.n = vn.xyz;
	//1 mean need to interpolate, compute from normals_tex
	else
	{
		vec3 vn0 = texelFetch(normals_tex, vn.x).xyz;
		vec3 vn1 = texelFetch(normals_tex, vn.y).xyz;
		vec3 vn2 = texelFetch(normals_tex, vn.z).xyz;

		rec.n = interpolate(vn0, vn1, vn2, rec.u, rec.v);
	}

	//Material mat;
	
	
	//rec.mat.tex_ind = texelFetch(mats_tex, 4 * mtl_ind + 3);
	rec.tex_ind = texelFetch(mats_tex, 4 * mtl_ind + 3);
	
	//vec2 texcoord = rec.texcoord;
	
	//rec.mat.emission = texelFetch(mats_tex, 4 * mtl_ind + 1);
	//rec.mat.specular = texelFetch(mats_tex, 4 * mtl_ind + 2);

	rec.emission = texelFetch(mats_tex, 4 * mtl_ind + 1);
	rec.specular = texelFetch(mats_tex, 4 * mtl_ind + 2);

	//if(rec.mat.tex_ind.x != -1)
	if(rec.tex_ind.x != -1)
	{
		ivec4 vt = texelFetch(triangles_tex, 3 * triangle_ind + 2);

		vec2 vt0 = texelFetch(texcoords_tex, vt.x).xy;
		vec2 vt1 = texelFetch(texcoords_tex, vt.y).xy;
		vec2 vt2 = texelFetch(texcoords_tex, vt.z).xy;

		rec.texcoord = interpolate(vt0, vt1, vt2, rec.u, rec.v);

		//rec.mat.albedo.xyz = pow(texture(albedo_textures, vec3(rec.texcoord, rec.mat.tex_ind.x)).xyz, vec3(2.2f));	
		rec.albedo.xyz = pow(texture(albedo_textures, vec3(rec.texcoord, rec.tex_ind.x)).xyz, vec3(2.2f));	
	}
	else
		//rec.mat.albedo = texelFetch(mats_tex, 4 * mtl_ind);
		rec.albedo = texelFetch(mats_tex, 4 * mtl_ind);

	//rec.mat = mat;
}

/*void compute_material(inout HitRecord rec)
{
	int mtl_ind = rec.mtl_ind;

	Material mat;

	mat.albedo = texelFetch(mats_tex, 4 * mtl_ind);
	mat.emission = texelFetch(mats_tex, 4 * mtl_ind + 1);
	mat.specular = texelFetch(mats_tex, 4 * mtl_ind + 2);
	mat.tex_ind = texelFetch(mats_tex, 4 * mtl_ind + 3);

	vec2 texcoord = rec.texcoord;
	
	if(mat.tex_ind.x != -1)
		mat.albedo.xyz = pow(texture(albedo_textures, vec3(texcoord, mat.tex_ind.x)).xyz, vec3(2.2f));

	rec.mat = mat;
}*/


bool trace_bvh(Ray r, out HitRecord rec)
{
	int stk[12];
	
	int ptr = 0;
	
	stk[ptr++] = -1;
	
	float t = inf;
	
	const vec3 invdir = 1.0f / (r.d);

	int ind = 0;
	
	
	while(ind > -1)
	{
		int n = ind;
		//vec3 info = texelFetch(bvh, n * 3 + 2).xyz;

		//int leftIndex = int(info.x);
		

		//float isLeaf = (info.z);

		vec4 bmin = texelFetch(bvh, 2 * n);
		vec4 bmax = texelFetch(bvh, 2 * n + 1);

		int leftIndex = int(bmin.w);

		if(bmax.w == 0)
		{
			//bfs
			int rightIndex = leftIndex + 1;
		

			//dfs
			//int rightIndex = int(info.y);

			float tl1, tl2;// th1;
			//float tl2;// th2;

			//vec3 o = r.o;

			float th1 = hit_bbox(r, texelFetch(bvh, leftIndex * 2).xyz, texelFetch(bvh, leftIndex * 2 + 1).xyz, invdir, tl1);
			float th2 = hit_bbox(r, texelFetch(bvh, rightIndex * 2).xyz, texelFetch(bvh, rightIndex * 2 + 1).xyz, invdir, tl2);


			//float th1 = hit_bbox(r, bmin.xyz, bmax.xyz, invdir, tl1);
			//float th2 = hit_bbox(r, bmin.xyz, bmax.xyz, invdir, tl2);

			bool left  = th1 > 0 && th1 >= tl1 && tl1 < t;//th1 >= 0
			bool right = th2 > 0 && th2 >= tl2 && tl2 < t;//th2 >= 0

			if(left)
			{
				ind = leftIndex;
				//int defered = rightIndex;

				if(right)
				{
					//if(tl1 > tl2)
					//{
					//	ind = rightIndex;
					//	defered = leftIndex;
					//}
					int node_offset = tl1 > tl2 ? 1 : 0;

					stk[ptr++] = ind + 1 - node_offset;//defered;

					ind += node_offset;
					//defered -= node_offset;
					

					
				}



				continue;
			}
			else if(right)
			{
				ind = rightIndex;
				continue;
			}


			/*if(left && right)
			{		
				//stk[++ind] = rightIndex;
				//stk[++ind] = leftIndex;

				//int defered = -1;
				
				ind = leftIndex;
				int defered = rightIndex;

				if (tl1 > tl2)
				{
					ind = rightIndex;
					defered = leftIndex;
				}
				//else
				//{
					//ind = leftIndex;
					//defered = rightIndex;
				//}
				stk[ptr++] = defered;
				continue;
			}			
			if(left)
			{
				ind = leftIndex;
				continue;
			}
			//them && right se bi mot soc den o giua
			//else if (rightHit >= 0.0 && right)//>= 0
			if(right)
			{
				ind = rightIndex;
				continue;
			}*/
			
		}	
		//if (isLeaf == 1)
		else
		{
			int range = int(bmax.w);//int(info.y);
			
			int end = leftIndex + range;
			for (int i = leftIndex; i < end; ++i) // Loop through indices
				hit_triangle(r, i, rec, t);

			//hit_triangle(r, leftIndex, rec, t);
			//		break;

			//hit_triangle_range(r, leftIndex, end, rec, t);
		}
		
		ind = stk[--ptr];
	}
	
	if(t < inf)
	{
		
		//vec4 triIndex = texelFetch(triangles_tex, 3 * rec.triangle_ind);
		rec.t = t;
		//rec.mtl_ind = int(triIndex.w);

		compute_normal_texcoord_and_material(rec);
		//compute_material(rec);

		return true;
	}
	return false;
}
//shared int stk[16];
bool hit_shadow(Ray r, float max_t)
{
	int stk[16];
	
	int ptr = 0;
	
	stk[ptr++] = -1;
	
	float t = inf;
	
	vec3 invdir = 1.0 / r.d;

	int ind = 0;

	while(ind > -1)
	{
		int n = ind;
		//vec3 info = texelFetch(bvh, n * 3 + 2).xyz;

		vec4 bmin = texelFetch(bvh, 2 * n);
		vec4 bmax = texelFetch(bvh, 2 * n + 1);

		int leftIndex = int(bmin.w);

		//int leftIndex = int(info.x);
		
		//float isLeaf = (info.z);

		//if (isLeaf == 1)
		if(bmax.w != 0)
		{
			int range = int(bmax.w);//int(info.y);

			
			for (int i = leftIndex; i < leftIndex + range; ++i) // Loop through indices
			//for (int i = leftIndex + range - 1; i >= leftIndex; --i)
				if(hit_triangle_no_rec(r, i, max_t))
					return true;
			//if(hit_triangle_no_rec(r, leftIndex, max_t))
			//	return true;		
		}
		else
		{
			//bfs
			int rightIndex = leftIndex + 1;

			//dfs
			//int rightIndex = int(info.y);

			//float tl1, tl2;
			//float leftHit = hit_bbox_2(r, texelFetch(bvh, leftIndex * 3 + 0).xyz, texelFetch(bvh, leftIndex * 3 + 1).xyz, invdir, tl1);
			//float rightHit = hit_bbox_2(r, texelFetch(bvh, rightIndex * 3 + 0).xyz, texelFetch(bvh, rightIndex * 3 + 1).xyz, invdir, tl2);
			
			//float leftHit = hit_bbox(r, texelFetch(bvh, leftIndex * 3 + 0).xyz, texelFetch(bvh, leftIndex * 3 + 1).xyz, invdir);
			//float rightHit = hit_bbox(r, texelFetch(bvh, rightIndex * 3 + 0).xyz, texelFetch(bvh, rightIndex * 3 + 1).xyz, invdir);
			//bool left = leftHit < t;
			//bool right = leftHit < t;

			float tl1, tl2;// th1;
			//float tl2;// th2;

			//vec3 o = r.o;
			//float th1 = hit_bbox(r, texelFetch(bvh, leftIndex * 3 + 0).xyz, texelFetch(bvh, leftIndex * 3 + 1).xyz, invdir, tl1);
			//float th2 = hit_bbox(r, texelFetch(bvh, rightIndex * 3 + 0).xyz, texelFetch(bvh, rightIndex * 3 + 1).xyz, invdir, tl2);

			
			float th1 = hit_bbox(r, texelFetch(bvh, leftIndex * 2).xyz, texelFetch(bvh, leftIndex * 2 + 1).xyz, invdir, tl1);
			float th2 = hit_bbox(r, texelFetch(bvh, rightIndex * 2).xyz, texelFetch(bvh, rightIndex * 2 + 1).xyz, invdir, tl2);

			//float th1 = hit_bbox(r, bmin.xyz, bmax.xyz, invdir, tl1);
			//float th2 = hit_bbox(r, bmax.xyz, bmin.xyz, invdir, tl2);

			bool left  = th1 >= 0 && th1 >= tl1 && tl1 <= max_t;
			bool right = th2 >= 0 && th2 >= tl2 && tl2 <= max_t;

			

			//bool left  = th1 >= 0 && th1 >= tl1 && tl1 <= max_t;
			//bool right = th2 >= 0 && th2 >= tl2 && tl2 <= max_t;

			//if (leftHit > 0.0 && rightHit > 0.0 && left && right)

			if(left)
			{
				ind = leftIndex;
				//int defered = rightIndex;

				if(right)
				{
					//if(tl1 > tl2)
					//{
					//	ind = rightIndex;
					//	defered = leftIndex;
					//}
					int node_offset = tl1 > tl2 ? 1 : 0;

					stk[ptr++] = ind + 1 - node_offset;//defered;

					ind += node_offset;
					//defered -= node_offset;
					

					
				}



				continue;
			}
			else if(right)
			{
				ind = rightIndex;
				continue;
			}
			/*if(left && right)
			{		
				//stk[++ind] = rightIndex;
				//stk[++ind] = leftIndex;

				int defered = -1;

				if (tl1 > tl2)
				{
					ind = rightIndex;
					defered = leftIndex;
				}
				else
				{
					ind = leftIndex;
					defered = rightIndex;
				}
				stk[ptr++] = defered;
				continue;
			}
			//else if (leftHit > 0.0 && left)//>= 0
			if(left)
			{
				ind = leftIndex;
				continue;
			}
			//else if (rightHit > 0.0 && right)//>= 0
			if(right)
			{
				ind = rightIndex;
				continue;
			}*/
		}	
		ind = stk[--ptr];
	}
	return false;
}

//for debug only
vec3 path_trace_albedo(inout Ray r)
{
	HitRecord rec;
	if(trace_bvh(r, rec))
	{
		
		int triangle_ind = rec.triangle_ind;

		//vec4 v = texelFetch(triangles_tex, 3 * triangle_ind);

		int mtl_ind = (rec.mtl_ind);//int(v.w);

		//vec3 albedo = texelFetch(mats_tex, 3 * mtl_ind).xyz;
		vec3 albedo = texelFetch(mats_tex, 4 * mtl_ind).xyz;

		return albedo;
	}
	return vec3(0.0f);
}

//sample triangle light position
vec3 sample_light_position(int li)
{
	vec3 p = texelFetch(lights_tex, 6 * li).xyz;
	vec3 u = texelFetch(lights_tex, 6 * li + 1).xyz;
	vec3 v = texelFetch(lights_tex, 6 * li + 2).xyz;

	float s = sqrt(rand());

	float b0 = 1.0f - s;
	float b1 = rand() * s;

	return p + b0 * u + b1 * v;
}

vec3 path_trace(inout Ray r)
{
	vec3 L = vec3(0.0f);
	vec3 T = vec3(1.0f);

	float prev_pdf = 1.0f;

	//Ray new_ray = r;
	bool is_specular = true;

	for(int i = 0; i < 3; ++i)
	{
		HitRecord rec;
		if(trace_bvh(r, rec))
		{
			float cos_incident = dot(r.d, rec.n);

			vec3 original_n = rec.n;

			if (cos_incident > 0)
					rec.n = -1.0f * rec.n;

			//int triangle_ind = rec.triangle_ind;

			//vec4 v = texelFetch(triangles_tex, 3 * triangle_ind);
			//ivec4 vt = texelFetch(triangles_tex, 3 * triangle_ind + 1);

			//int mtl_ind = rec.mtl_ind;//int(v.w);

			//vec4 albedo = texelFetch(mats_tex, 3 * mtl_ind);
			//vec4 emission = texelFetch(mats_tex, 3 * mtl_ind + 1);

			//vec4 albedo = rec.mat.albedo;
			//vec4 emission = rec.mat.emission;

			vec4 emission = rec.emission;

			if(emission.w != -1)
			{
				if(is_specular)
				{
					L += T * emission.xyz;
					return L;
				}
				else
				{
					vec3 light_direction = r.d * rec.t;

					float length = length(light_direction);

					light_direction = normalize(light_direction);

					float cos_light = (-1.0f * dot(light_direction, rec.n));//ko can abs

					float length2 = length * length;

					int light_index = int(emission.w);//int(vt.w);

					vec3 light_area_pdf = texelFetch(lights_tex, 6 * light_index + 5).xyz;

					//area_pdf.x = area of this light source
					//area_pdf.y = pdf of this light source

					float pdf_light = length2 / (light_area_pdf.x * cos_light) * light_area_pdf.y;

					float mis_weight = power_heuristic(prev_pdf, pdf_light);

					L += T * emission.xyz * mis_weight;

					return L;
				}
			}

			vec3 hit_point = r.o + r.d * rec.t + rec.n * 0.0002f;

			//vec4 specular = rec.mat.specular;
			vec4 specular = rec.specular;

			//vec4 specular = texelFetch(mats_tex, 3 * mtl_ind + 2);

			//not specular material
			if(specular.w == 0)
			{
				int light_index = int(rand() * numLights);
				
				vec3 light_position = sample_light_position(light_index);

				vec3 light_direction = light_position - hit_point;

				float length = length(light_direction);

				float ilength = 1.0f / length;

				light_direction *= ilength;

				float cos_mtl = dot(light_direction, original_n);

				//float cos_light = dot(light_direction, rec.n);

				Ray light_ray = Ray(hit_point, light_direction);

				vec3 light_normal = texelFetch(lights_tex, 6 * light_index + 3).xyz;


				//neu de cos_light nay se gay ra nhieu vung firefly
				//float cos_light = abs((dot(light_direction, light_normal)));

				//float cos_light = (-1.0f * (dot(light_direction, light_normal)));

				float cos_light = ((dot(light_direction, light_normal)));

				if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_shadow(light_ray, length - eps))
				//if(!hit_shadow(light_ray, length - eps))
				{					
					vec3 light_emission = texelFetch(lights_tex, 6 * light_index + 4).xyz;
					
					

					vec3 light_area_pdf = texelFetch(lights_tex, 6 * light_index + 5).xyz;

					//de cos_light o day gay ra nhieu vung sang ao
					//float cos_light = abs(-1.0f * dot(light_direction, light_normal));

					

					

					//float length2 = length * length;

					float pdf_light = (length * length) / (light_area_pdf.x * -cos_light) * light_area_pdf.y;

					vec3 bsdf_eval = diffuse_bsdf(r, rec, light_direction);

					float bsdf_pdf = diffuse_pdf(r, rec, light_direction);

					float mis_weight = power_heuristic(pdf_light, bsdf_pdf);

					//L += T * 2 * bsdf_eval * mis_weight / pdf_light;
					
					//L += T * 3 * bsdf_eval * mis_weight / pdf_light;

					L += T * light_emission * bsdf_eval * mis_weight / pdf_light;
					//float cos_mtl = dot(light_direction, rec.n);
					
				}
			}

			//Sample for new ray
			vec3 sample_direction = diffuse_sample(r, rec.n);
			vec3 new_hit_point = r.o + r.d * rec.t + rec.n * 0.0002f;

			vec3 sample_eval = diffuse_bsdf(r, rec, sample_direction);

			float bsdf_pdf = diffuse_pdf(r, rec, sample_direction);

			T *= sample_eval;

			prev_pdf = bsdf_pdf;

			is_specular = false;

			r = Ray(new_hit_point, sample_direction);
		}
		else
			return L;
	}
	return L;
}

void main()
{
	seed = gl_FragCoord.xy;

	float r1 = 2.0f * rand();
	float r2 = 2.0f * rand();

	vec2 jitter;

	jitter.x = r1 < 1.0 ? sqrt(r1) - 1.0f : 1.0f - sqrt(2.0f - r1);
	jitter.y = r2 < 1.0 ? sqrt(r2) - 1.0f : 1.0f - sqrt(2.0f - r2);
	jitter /= (screenResolution * 0.5f);

	vec2 d = (2.0f * tex - 1.0f) + jitter;

	float tan_fov = tan(camera.fov * 0.5f);
	d.x *= screenResolution.x / screenResolution.y * tan_fov;
	d.y *= tan_fov;

	vec3 rayDir = normalize(d.x * camera.right + d.y * camera.up + camera.forward);

	Ray ray = Ray(camera.position, rayDir);

	//vec3 pixelColor = Path_Tracer2(ray);
	//vec3 pixelColor = Path_Tracer(ray);
	//vec3 pixelColor = brute_force(ray);//Path_Tracer(ray);

	//vec3 pixelColor = path_trace_albedo(ray);

	vec3 accumulate_color = texture(accumulate_tex, tex).xyz;

	vec3 pixelColor = path_trace(ray);

	color = pixelColor + accumulate_color;
}


