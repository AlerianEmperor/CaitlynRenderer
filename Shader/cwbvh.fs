R(#version 460

in vec2 tex;

out vec3 color;
uniform sampler2D accumulate_tex;
uniform samplerBuffer vertices_tex;
uniform samplerBuffer normals_tex;
uniform samplerBuffer texcoords_tex;
uniform isamplerBuffer triangles_tex; //isamplerBuffer: reduce float to int conversion, result in 2.8% speed up!
uniform samplerBuffer mats_tex; 
uniform samplerBuffer lights_tex;
//uniform samplerBuffer bvh;
uniform samplerBuffer bvh8;
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
	bmin = (bmin - r.o) * invdir;
	bmax = (bmax - r.o) * invdir;
		
	vec3 tmax = max(bmax, bmin);
	bmin = min(bmax, bmin);

	float th = min(tmax.x, min(tmax.y, tmax.z));
	
	tl = max(bmin.x, max(bmin.y, bmin.z));

	return th;
}



float power_heuristic(float a, float b)
{
	float t = a * a;
	return t / (b * b + t);
}


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
	
	int mtl_ind;
	//int mtl_ind;
};

vec3 cosine_hemisphere_sampling()
{
	float u1 = rand();
	float u2 = rand();
	
	float r = sqrt(u1);
	float phi = pi2 * u2;


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

//CWBVH
//get oct_inv4 of ray direction
uint get_oct_inv4(in vec3 d)
{
	return (d.x < 0.0f ? 0 : 0x04040404) |
		   (d.y < 0.0f ? 0 : 0x02020202) |
		   (d.z < 0.0f ? 0 : 0x01010101);
}

struct BVH8node
{
	vec4 node0;
	vec4 node1;
	vec4 node2;
	vec4 node3;
	vec4 node4;
};

uint extract_byte(uint x, uint i)
{
	return (x >> (8 * i)) & 0xff;
}

uint sign_extend_s8x4(uint x)
{
	return ((x >> 7) & 0x01010101) * 0xff;
}

#define LOCAL_STACK_SIZE 16

uint bvh8_node_intersect(in Ray r, vec3 inv_direction, uint oct_inv4, float max_t, BVH8node node)
{
	vec3 p = node.node0.xyz;

	uint e_imask = floatBitsToUint(node.node0.w);
	uint e_x = extract_byte(e_imask, 0);
	uint e_y = extract_byte(e_imask, 1);
	uint e_z = extract_byte(e_imask, 2);

	vec3 adjust_ray_direction_inv = vec3(uintBitsToFloat(e_x << 23) * inv_direction.x,
										 uintBitsToFloat(e_y << 23) * inv_direction.y,
										 uintBitsToFloat(e_z << 23) * inv_direction.z);
	
	vec3 adjust_ray_origin = (p - r.o) * inv_direction;

	uint hit_mask = 0;

	for(int i = 0; i < 2; ++i)
	{
		uint meta4 = floatBitsToUint(i == 0 ? node.node1.z : node.node1.w);

		uint is_inner4 = (meta4 & (meta4 << 1)) & 0x10101010;
		uint inner_mask4 = sign_extend_s8x4(is_inner4 << 3);
		
		uint bit_index4 = (meta4 ^ (oct_inv4 & inner_mask4)) & 0x1F1F1F1F;

		uint child_bits4 = (meta4 >> 5) & 0x07070707;

		//near far
		uint q_lo_x = floatBitsToUint(i == 0 ? node.node2.x : node.node2.y);
		uint q_hi_x = floatBitsToUint(i == 0 ? node.node2.z : node.node2.w);

		uint q_lo_y = floatBitsToUint(i == 0 ? node.node3.x : node.node3.y);
		uint q_hi_y = floatBitsToUint(i == 0 ? node.node3.z : node.node3.w);

		uint q_lo_z = floatBitsToUint(i == 0 ? node.node4.x : node.node4.y);
		uint q_hi_z = floatBitsToUint(i == 0 ? node.node4.z : node.node4.w);

		uint x_min = r.d.x < 0.0f ? q_hi_x : q_lo_x;
		uint x_max = r.d.x < 0.0f ? q_lo_x : q_hi_x;

		uint y_min = r.d.y < 0.0f ? q_hi_y : q_lo_y;
		uint y_max = r.d.y < 0.0f ? q_lo_y : q_hi_y;

		uint z_min = r.d.z < 0.0f ? q_hi_z : q_lo_z;
		uint z_max = r.d.z < 0.0f ? q_lo_z : q_hi_z;

		//intersect 4 children
		//1 node have 2 children
		//each children have another 4 children or childs! (childs sound like child ass! that's why it's called children!)
		for(int j = 0; j < 4; ++j)
		{
			vec3 tmin3 = vec3(float(extract_byte(x_min, j)), float(extract_byte(y_min, j)),float(extract_byte(z_min, j)));
			vec3 tmax3 = vec3(float(extract_byte(x_max, j)), float(extract_byte(y_max, j)),float(extract_byte(z_max, j)));

			tmin3 = tmin3 * adjust_ray_direction_inv + adjust_ray_origin;
			tmax3 = tmax3 * adjust_ray_direction_inv + adjust_ray_origin;

			float tmin = max(max(tmin3.x, tmin3.y), tmin3.z);
			float tmax = max(max(tmax3.x, tmax3.y), tmax3.z);

			if(tmax >= 0.0f && tmax < max_t && tmin <= tmax)
			{
				uint child_bits = extract_byte(child_bits4, j);
				uint bit_index = extract_byte(bit_index4, j);
				hit_mask |= child_bits << bit_index;
			}
		}
	}
	return hit_mask;
}

bool hit_bvh8(in Ray r, inout HitRecord rec)
{
	float max_t = inf;

	uvec2 local_stack[LOCAL_STACK_SIZE];

	int stack_size = 0;

	uint oct_inv4 = get_oct_inv4(r.d);

	uvec2 current_group = uvec2(0, 0x80000000);

	vec3 inv_direction = 1.0f / r.d;

	while(stack_size > 0 || current_group.y != 0)
	{
		uvec2 triangle_group;

		if((current_group.y & 0xff000000) != 0)
		{
			uint hits_imask = current_group.y;
			int child_index_offset = findMSB(hits_imask);
			uint child_index_base = current_group.x;

			current_group.y &= ~(1 << child_index_offset);

			if((current_group.y & 0xff000000) != 0)
				local_stack[stack_size++] = current_group;

			uint slot_index = (child_index_offset - 24) ^ (oct_inv4 & 0xff);
			uint relative_index = bitCount(hits_imask & ~(0xffffffff << slot_index));

			uint child_node_index = child_index_base + relative_index;

			BVH8node node;

			node.node0 = texelFetch(bvh8, int(child_node_index * 5));
			node.node1 = texelFetch(bvh8, int(child_node_index * 5 + 1));
			node.node2 = texelFetch(bvh8, int(child_node_index * 5 + 2));
			node.node3 = texelFetch(bvh8, int(child_node_index * 5 + 3));
			node.node4 = texelFetch(bvh8, int(child_node_index * 5 + 4));

			//uint bvh8_node_intersect(in Ray r, vec3 inv_direction, uint oct_inv4, float max_t, BVH8node node)

			uint hitmask = bvh8_node_intersect(r, inv_direction, oct_inv4, max_t, node);

			uint imask = extract_byte(floatBitsToUint(node.node0.w), 3);

			current_group.x = floatBitsToUint(node.node1.x);
			triangle_group.x = floatBitsToUint(node.node1.y);

			current_group.y = (hitmask & 0xff000000) | imask;
			triangle_group.y = (hitmask & 0x00ffffff);
		}
		else
		{
			triangle_group = current_group;
			current_group = uvec2(0);
		}

		while(triangle_group.y != 0)
		{
			int triangle_ind = findMSB(triangle_group.y);
			triangle_group.y &= ~(1 << triangle_ind);

			int tri_idx = int(triangle_group.x + triangle_ind);

			//bool hit_triangle(Ray r, int i, inout HitRecord rec, inout float t)

			hit_triangle(r, tri_idx, rec, max_t);
		}

		if((current_group.y & 0xff000000) == 0)
		{
			if(stack_size == 0)
				break;
			current_group = local_stack[stack_size--];
		}
	}
	if(max_t < inf)
	{
		rec.t = max_t;
		
		compute_normal_texcoord_and_material(rec);
	
		return true;
	}
	return false;
}

bool hit_bvh8_shadow(in Ray r, float max_t)
{
	uvec2 local_stack[LOCAL_STACK_SIZE];

	int stack_size = 0;

	uint oct_inv4 = get_oct_inv4(r.d);

	uvec2 current_group = uvec2(0, 0x80000000);

	vec3 inv_direction = 1.0f / r.d;

	while(stack_size > 0 || current_group.y != 0)
	{
		uvec2 triangle_group;

		if((current_group.y & 0xff000000) != 0)
		{
			uint hits_imask = current_group.y;
			int child_index_offset = findMSB(hits_imask);
			uint child_index_base = current_group.x;

			current_group.y &= ~(1 << child_index_offset);

			if((current_group.y & 0xff000000) != 0)
				local_stack[stack_size++] = current_group;

			uint slot_index = (child_index_offset - 24) ^ (oct_inv4 & 0xff);
			uint relative_index = bitCount(hits_imask & ~(0xffffffff << slot_index));

			uint child_node_index = child_index_base + relative_index;

			BVH8node node;

			node.node0 = texelFetch(bvh8, int(child_node_index * 5));
			node.node1 = texelFetch(bvh8, int(child_node_index * 5 + 1));
			node.node2 = texelFetch(bvh8, int(child_node_index * 5 + 2));
			node.node3 = texelFetch(bvh8, int(child_node_index * 5 + 3));
			node.node4 = texelFetch(bvh8, int(child_node_index * 5 + 4));

			//uint bvh8_node_intersect(in Ray r, vec3 inv_direction, uint oct_inv4, float max_t, BVH8node node)

			uint hitmask = bvh8_node_intersect(r, inv_direction, oct_inv4, max_t, node);

			uint imask = extract_byte(floatBitsToUint(node.node0.w), 3);

			current_group.x = floatBitsToUint(node.node1.x);
			triangle_group.x = floatBitsToUint(node.node1.y);

			current_group.y = (hitmask & 0xff000000) | imask;
			triangle_group.y = (hitmask & 0x00ffffff);
		}
		else
		{
			triangle_group = current_group;
			current_group = uvec2(0);
		}

		while(triangle_group.y != 0)
		{
			int triangle_ind = findMSB(triangle_group.y);
			triangle_group.y &= ~(1 << triangle_ind);

			int tri_idx = int(triangle_group.x + triangle_ind);

			if(hit_triangle_no_rec(r, tri_idx, max_t))
				return true;
		}

		if((current_group.y & 0xff000000) == 0)
		{
			if(stack_size == 0)
				break;
			current_group = local_stack[stack_size--];
		}
	}

	return false;
}


/*
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
			
		}	
		ind = stk[--ptr];
	}
	return false;
}
*/

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
		//if(trace_bvh(r, rec))
		if(hit_bvh8(r, rec))
		{
			float cos_incident = dot(r.d, rec.n);

			vec3 original_n = rec.n;

			if (cos_incident > 0)
					rec.n = -1.0f * rec.n;

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

				//if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_shadow(light_ray, length - eps))
				//if(!hit_shadow(light_ray, length - eps))
				if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_bvh8_shadow(light_ray, length - eps))
				{					
					vec3 light_emission = texelFetch(lights_tex, 6 * light_index + 4).xyz;
					
					

					vec3 light_area_pdf = texelFetch(lights_tex, 6 * light_index + 5).xyz;

					float pdf_light = (length * length) / (light_area_pdf.x * -cos_light) * light_area_pdf.y;

					vec3 bsdf_eval = diffuse_bsdf(r, rec, light_direction);

					float bsdf_pdf = diffuse_pdf(r, rec, light_direction);

					float mis_weight = power_heuristic(pdf_light, bsdf_pdf);
					
					L += T * light_emission * bsdf_eval * mis_weight / pdf_light;
										
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

	

	vec3 accumulate_color = texture(accumulate_tex, tex).xyz;

	vec3 pixelColor = path_trace(ray);

	color = pixelColor + accumulate_color;
})

