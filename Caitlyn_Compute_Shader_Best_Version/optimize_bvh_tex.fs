#version 460

in vec2 tex;

out vec3 color;
uniform sampler2D accumulate_tex;
//uniform samplerBuffer vertices_tex;

uniform int tree_height;
int bvh_height = tree_height;


layout(std430, binding = 1) restrict readonly buffer vertices_buf
{
	vec4 vertices[];
};

//uniform samplerBuffer normals_tex;

layout(std430, binding = 2) buffer normals_buf
{
	vec4 normals[];
};

//uniform samplerBuffer texcoords_tex;

layout(std430, binding = 3) buffer texcoords_buf
{
	vec2 texcoords[];
};

//uniform isamplerBuffer triangles_tex; //isamplerBuffer: reduce float to int conversion, result in 2.8% speed up!
layout(std430, binding = 4) restrict readonly buffer triangles_buf
{
	ivec4 triangles[];
};

//uniform samplerBuffer mats_tex; 

layout(std430, binding = 5) buffer mats_buf
{
	vec4 materials[];
};

struct Light4
{
    vec4 p;
    vec4 u;
    vec4 v;
    vec4 e;
    vec4 n;
    vec4 area_pdf;
	//vec4 pad; vec4 pad2;
};

//uniform samplerBuffer lights_tex;
layout(std430, binding = 6) buffer light_buf
{
	Light4 lights[];
};
//uniform samplerBuffer bvh;

layout(std430, binding = 7) readonly buffer bvh_buf
{
	vec4 bvh[];
};
uniform sampler2DArray albedo_textures;

#define pi 3.1415926
#define pi2 6.2831853
#define ipi 1.0f / pi


#define VERTICES vertices_buf_block.vertices
#define NORMALS normals_buf_block.normals
#define TEXCOORDS texcoords_buf_block.texcoords
#define TRIANGLES triangles_buf_block.triangles
#define MATERIALS mats_buf_block.materials
#define BVH bvh_buf_block.bvh
#define LIGHTS light_buf_block.lights

struct Ray { vec3 o; vec3 d; };
struct Camera { vec3 up; vec3 right; vec3 forward; vec3 position; float fov; float focalDist; float aperture; };


//struct Light { vec3 position; vec3 emission; vec3 u; vec3 v; vec3 area; };


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

	uint state = floatBitsToUint(seed.x * 1234567.0 + seed.y * 7654321.0);
	
	return fract(sin(dot(seed, vec2(12.9898, 78.233))) * 43758.5453);
}

void onb(in vec3 n, out vec3 b1, out vec3 b2)
{
    if (n.z < -0.999f) {
        b1 = vec3(0, -1, 0);
        b2 = vec3(-1, 0, 0);
    } else {
        float a = 1.0/(1.0+n.z);
        float b = -n.x*n.y*a;
        b1 = vec3(1.0+b, b, -n.x);
        b2 = vec3(b, 1.0+b, -n.y);
    }
}

float hit_bbox(Ray r, vec3 bmin, vec3 bmax, vec3 invdir, out float tl)
{
	/*// shift origin once
    vec3 t1 = (bmin - r.o) * invdir;
    vec3 t2 = (bmax - r.o) * invdir;

    // get min / max along each axis
    float tx1 = min(t1.x, t2.x), tx2 = max(t1.x, t2.x);

    float ty1 = min(t1.y, t2.y), ty2 = max(t1.y, t2.y);

    float tz1 = min(t1.z, t2.z), tz2 = max(t1.z, t2.z);

    // merge intervals
    tl = max(tx1, max(ty1, tz1));
    return min(tx2, min(ty2, tz2));*/

	// 1. Tính toán song song trên c? 3 tr?c cùng lúc b?ng Vector
    vec3 t1 = (bmin - r.o) * invdir;
    vec3 t2 = (bmax - r.o) * invdir;

    // 2. Gom c?m c?c ti?u và c?c ??i c?a t?ng c?p m?t ph?ng (Không dùng if-else ?n)
    // GPU th?c hi?n hàm min/max c?a vec3 b?ng ?úng 2 chu k? l?nh ph?n c?ng (SIMD)
    vec3 tmin = min(t1, t2);
    vec3 tmax = max(t1, t2);

    // 3. Tìm kho?ng giao nhau l?n nh?t c?a l?i vào (tl) và nh? nh?t c?a l?i ra (th)
    // Dùng c?u trúc max/min l?ng nhau ki?u vector giúp trình biên d?ch t?i ?u hóa c?c t?t
    tl = max(tmin.x, max(tmin.y, tmin.z));
    float th = min(tmax.x, min(tmax.y, tmax.z));

    // Tr? v? kho?ng cách l?i ra. ?? ki?m tra tia có ?âm trúng h?p không, 
    // ?i?u ki?n ? hàm trace_bvh c?a b?n ch? c?n là: (tl <= th) && (th > 0.0)
    return th; 
}


float power_heuristic(float a, float b)
{
	float t = a * a;
	return t / (b * b + t);
}

float power_heuristic_inv_parameter(float inv_a2, float b)
{
	return 1.0f / (b * b * inv_a2 + 1.0f);
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
	//ivec4 triIndex = texelFetch(triangles_tex, 3 * i);
	
	ivec4 triIndex = triangles[3 * i];

	vec3 v0 = vertices[triIndex.x].xyz; // int(triIndex.x)
	vec3 v1 = vertices[triIndex.y].xyz - v0;
	vec3 v2 = vertices[triIndex.z].xyz - v0;
	
	//v1 -= v0;
	//v2 -= v0;

	//vec3 pv = cross(r.d, v2);
	//vec3 qv = cross(tv, v1);

	vec3 pv = vec3(
    r.d.y * v2.z - r.d.z * v2.y,
    r.d.z * v2.x - r.d.x * v2.z,
    r.d.x * v2.y - r.d.y * v2.x);

	vec3 tv = r.o - v0;
	
	vec3 qv = vec3(
    tv.y * v1.z - tv.z * v1.y,
    tv.z * v1.x - tv.x * v1.z,
    tv.x * v1.y - tv.y * v1.x
);

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


bool hitTriangleNoRec(in Ray r, int i, inout float t)
{
    ivec4 triIndex = triangles[3 * i];

    // Load vertices
    vec3 v0 = vertices[triIndex.x].xyz;
    vec3 v1 = vertices[triIndex.y].xyz - v0;
    vec3 v2 = vertices[triIndex.z].xyz - v0;

    // Convert to edges
    //v1 -= v0;
    //v2 -= v0;

    // Cross( r.d, v2 ) ? pv
    vec3 pv = vec3(
        r.d.y * v2.z - r.d.z * v2.y,
        r.d.z * v2.x - r.d.x * v2.z,
        r.d.x * v2.y - r.d.y * v2.x
    );

    float det = v1.x * pv.x + v1.y * pv.y + v1.z * pv.z;
    if (det == 0.0) return false;

    float inv_det = 1.0 / det;

    // tv = r.o - v0
    vec3 tv = r.o - v0;

    // u = dot(tv, pv) * inv_det
    float u = (tv.x * pv.x + tv.y * pv.y + tv.z * pv.z) * inv_det;
    if (u < 0.0 || u > 1.0) return false;

    // qv = cross(tv, v1)
    vec3 qv = vec3(
        tv.y * v1.z - tv.z * v1.y,
        tv.z * v1.x - tv.x * v1.z,
        tv.x * v1.y - tv.y * v1.x
    );

    // v = dot(r.d, qv) * inv_det
    float v = (r.d.x * qv.x + r.d.y * qv.y + r.d.z * qv.z) * inv_det;
    if (v < 0.0 || u + v > 1.0) return false;

    // t = dot(v2, qv) * inv_det
    float tHit = (v2.x * qv.x + v2.y * qv.y + v2.z * qv.z) * inv_det;
    if (tHit <= 0.0 || tHit >= t) return false;

    // Commit hit
    //t = tHit;
    //rec.u = u;
   //rec.v = v;
    //rec.triangle_ind = i;
   // rec.mtl_ind = triIndex.w;

    return true;
}

bool hit_triangle_no_rec(Ray r, int i, float max_t)
{
	//ivec4 triIndex = texelFetch(triangles_tex, 3 * i);
	
	ivec4 triIndex = triangles[3 * i];

	vec3 v0 = vertices[triIndex.x].xyz; // int(triIndex.x)
	vec3 v1 = vertices[triIndex.y].xyz;
	vec3 v2 = vertices[triIndex.z].xyz;


	//vec3 v0 = texelFetch(vertices_tex, triIndex.x).xyz; // int(triIndex.x)
	//vec3 v1 = texelFetch(vertices_tex, triIndex.y).xyz;
	//vec3 v2 = texelFetch(vertices_tex, triIndex.z).xyz;

	//vec3 e0 = v1 - v0;
	//vec3 e1 = v2 - v0;
	v1 -= v0;
	v2 -= v0;

	/*vec3 pv = cross(r.d, v2);
	//float det = dot(e0, pv);

	vec3 tv = r.o - v0;
	vec3 qv = cross(tv, v1);*/

	vec3 pv = vec3(
    r.d.y * v2.z - r.d.z * v2.y,
    r.d.z * v2.x - r.d.x * v2.z,
    r.d.x * v2.y - r.d.y * v2.x);

	vec3 tv = r.o - v0;
	
	//vec3 qv = cross(tv, v1);

	vec3 qv = vec3(
    tv.y * v1.z - tv.z * v1.y,
    tv.z * v1.x - tv.x * v1.z,
    tv.x * v1.y - tv.y * v1.x);

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

	
	//ivec4 vn = texelFetch(triangles_tex, 3 * triangle_ind + 1);

	ivec4 vn = triangles[3 * triangle_ind + 1];
	
	//0 mean no need to interpolate, just use directly
	if(vn.w == 0)
		rec.n = vn.xyz;
	//1 mean need to interpolate, compute from normals_tex
	else
	{
		//vec3 vn0 = texelFetch(normals_tex, vn.x).xyz;
		//vec3 vn1 = texelFetch(normals_tex, vn.y).xyz;
		//vec3 vn2 = texelFetch(normals_tex, vn.z).xyz;

		vec3 vn0 = normals[vn.x].xyz;
		vec3 vn1 = normals[vn.y].xyz;
		vec3 vn2 = normals[vn.z].xyz;

		rec.n = interpolate(vn0, vn1, vn2, rec.u, rec.v);
	}
	
	rec.tex_ind = materials[4 * mtl_ind + 3];
	
	rec.emission = materials[4 * mtl_ind + 1];
	rec.specular = materials[4 * mtl_ind + 2];

	//if(rec.mat.tex_ind.x != -1)
	if(rec.tex_ind.x != -1)
	{
		//ivec4 vt = texelFetch(triangles_tex, 3 * triangle_ind + 2);

		ivec4 vt = triangles[3 * triangle_ind + 2];

		//vec2 vt0 = texelFetch(texcoords_tex, vt.x).xy;
		//vec2 vt1 = texelFetch(texcoords_tex, vt.y).xy;
		//vec2 vt2 = texelFetch(texcoords_tex, vt.z).xy;

		vec2 vt0 = texcoords[vt.x].xy;//texelFetch(texcoords_tex, vt.x).xy;
		vec2 vt1 = texcoords[vt.y].xy;//texelFetch(texcoords_tex, vt.y).xy;
		vec2 vt2 = texcoords[vt.z].xy;//texelFetch(texcoords_tex, vt.z).xy;

		rec.texcoord = interpolate(vt0, vt1, vt2, rec.u, rec.v);

		//rec.mat.albedo.xyz = pow(texture(albedo_textures, vec3(rec.texcoord, rec.mat.tex_ind.x)).xyz, vec3(2.2f));	
		//rec.albedo.xyz = pow(texture(albedo_textures, vec3(rec.texcoord, rec.tex_ind.x)).xyz, vec3(2.2f));	
		rec.albedo.xyz = texture(albedo_textures, vec3(rec.texcoord, rec.tex_ind.x)).xyz;
	}
	else
		//rec.mat.albedo = texelFetch(mats_tex, 4 * mtl_ind);
		//rec.albedo = texelFetch(mats_tex, 4 * mtl_ind);
		rec.albedo = materials[4 * mtl_ind];

	//rec.mat = mat;
}

bool trace_bvh_dfs_final_2(Ray r, out HitRecord rec)
{
    int stk[32]; // Stack lưu Node ID (0, 1, 2...)
    int ptr = 0;
    
    float t = inf;
    const vec3 invdir = 1.0f / (r.d);
    int ind = 0; // Bắt đầu từ Root ID = 0
    
    while (ind != -1)
    {
        // Mỗi Node chiếm 2 phần tử vec4 liên tiếp: bvh[ind * 2] và bvh[ind * 2 + 1]
        int n = ind << 1; 
        vec4 bmin = bvh[n];
        vec4 bmax = bvh[n + 1];

        int range = int(bmax.w);

        if (range == 0) // --- NÚT TRONG ---
        {
            int leftIndex  = ind + 1;     // Theo DFS, con trái luôn là nút kế tiếp
            int rightIndex = int(bmin.w); // Con phải lấy từ chỉ số nhảy được CPU tính sẵn

            float tl1, th1, tl2, th2;
            th1 = hit_bbox(r, bvh[leftIndex * 2].xyz,  bvh[leftIndex * 2 + 1].xyz,  invdir, tl1);
            th2 = hit_bbox(r, bvh[rightIndex * 2].xyz, bvh[rightIndex * 2 + 1].xyz, invdir, tl2);

            bool left  = (tl1 <= th1) && (th1 > 0.0) && (tl1 <= t);
            bool right = (tl2 <= th2) && (th2 > 0.0) && (tl2 <= t);

            if (left && right)
            {
                // Ray Sorting: Ưu tiên đi vào hộp gần hơn trước để tối ưu cắt giảm nhánh
                if (tl1 < tl2)
                {
                    stk[ptr++] = rightIndex; // Hoãn con phải
                    ind = leftIndex;         // Đi xuống con trái
                }
                else
                {
                    stk[ptr++] = leftIndex;  // Hoãn con trái
                    ind = rightIndex;        // Đi xuống con phải
                }
                continue;
            }
            else if (left)
            {
                ind = leftIndex;
                continue;
            }
            else if (right)
            {
                ind = rightIndex;
                continue;
            }
        }   
        else // --- NÚT LÁ ---
        {
            int start = int(bmin.w);
            int end   = start + range;

            for (int i = start; i < end; ++i) 
            {
                hit_triangle(r, i, rec, t); 
            }
        }
        
        // Quay lui (Pop Stack)
        ind = (ptr > 0) ? stk[--ptr] : -1;
    }
    
    if (t < inf)
    {
        rec.t = t;
        compute_normal_texcoord_and_material(rec);
        return true;
    }
    return false;
}

bool hit_shadow_dfs_final(Ray r, float max_t)
{
    // Kích thước stack 32 là thoải mái cho DFS trên GPU
    int stk[32];
    int ptr = 0;
    
    const vec3 invdir = 1.0f / r.d;

    // Bắt đầu từ Node Root (vị trí phần tử 0 trong mảng bvh)
    int ind = 0;

    while (ind != -1)
    {
        // Đọc dữ liệu của Node hiện tại
        int n = ind << 1; 
        vec4 bmin = bvh[n];
        vec4 bmax = bvh[n + 1];

        // w của bmax lưu range: 0 là Node trong, >0 là Node lá
        int range = int(bmax.w);

        if (range == 0) // --- NÚT TRONG (INTERNAL NODE) ---
        {
            // Theo cấu trúc DFS: Nút con trái luôn ở ngay kế tiếp nút cha
            int leftNodeArrayIdx  = ind + 1;//ind + 2; 
            
            // Nút con phải được tính từ chỉ số Node ID lưu trong bmin.w
            int rightNodeArrayIdx = int(bmin.w); 

            float tl1, th1, tl2, th2;
            
            // Tính toán giao cắt BBox cho cả 2 con
            th1 = hit_bbox(r, bvh[leftNodeArrayIdx * 2].xyz,  bvh[leftNodeArrayIdx * 2 + 1].xyz,  invdir, tl1);
            th2 = hit_bbox(r, bvh[rightNodeArrayIdx * 2].xyz, bvh[rightNodeArrayIdx * 2 + 1].xyz, invdir, tl2);

            // Điều kiện đâm trúng hộp và nằm trong khoảng max_t của tia shadow
            bool left  = (tl1 <= th1) && (th1 > 0.0) && (tl1 <= max_t);
            bool right = (tl2 <= th2) && (th2 > 0.0) && (tl2 <= max_t);

            if (left && right)
            {
                // Đối với Shadow Ray: Không cần quan tâm con nào gần hơn.
                // Cứ đẩy con phải vào stack và ưu tiên đi tiếp xuống con trái luôn để giảm chỉ thị rẽ nhánh.
                stk[ptr++] = rightNodeArrayIdx;
                ind = leftNodeArrayIdx;
                continue;
            }
            else if (left)
            {
                ind = leftNodeArrayIdx;
                continue;
            }
            else if (right)
            {
                ind = rightNodeArrayIdx;
                continue;
            }
        }
        else // --- NÚT LÁ (LEAF NODE) ---
        {
            int start = int(bmin.w);
            int end   = start + range;

            for (int i = start; i < end; ++i) 
            {
                // Nếu trúng bất kỳ tam giác nào, kết luận có bóng đổ và thoát luôn lập tức (Early Out)
                if (hitTriangleNoRec(r, i, max_t))
                    return true;
            }
        }	
        
        // POP STACK: Lấy nút hoãn ra duyệt tiếp, nếu stack trống thì gán -1 để thoát vòng lặp
        ind = (ptr > 0) ? stk[--ptr] : -1;
    }
    
    // Đã duyệt hết các nhánh tiềm năng mà không bị chắn => Không có bóng đổ
    return false;
}
//sample triangle light position
vec3 sample_light_position(int li)
{
	//vec3 p = texelFetch(lights_tex, 6 * li).xyz;
	//vec3 u = texelFetch(lights_tex, 6 * li + 1).xyz;
	//vec3 v = texelFetch(lights_tex, 6 * li + 2).xyz;

	Light4 l = lights[li];

	vec3 p = l.p.xyz;//lights[6 * li].xyz;
	vec3 u = l.u.xyz;//lights[6 * li + 1].xyz;
	vec3 v = l.v.xyz;//lights[6 * li + 2].xyz;

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
		//if(trace_bvh_final(r, rec))
		//if(trace_bvh(r, rec))
		//if(trace_bvh_dfs_final(r, rec))
		if(trace_bvh_dfs_final_2(r, rec))
		{
			float cos_incident = dot(r.d, rec.n);

			vec3 original_n = rec.n;

			if (cos_incident > 0)
					rec.n = -rec.n;

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

					float cos_light = (-dot(light_direction, rec.n));//ko can abs

					float length2 = length * length;

					int light_index = int(emission.w);//int(vt.w);

					//vec3 light_area_pdf = texelFetch(lights_tex, 6 * light_index + 5).xyz;

					Light4 l = lights[light_index];
					vec3 light_area_pdf = l.area_pdf.xyz;
					//vec3 light_area_pdf = lights[6 * light_index + 5].xyz;

					//area_pdf.x = area of this light source
					//area_pdf.y = pdf of this light source

					float pdf_light = length2 / (light_area_pdf.x * cos_light) * light_area_pdf.y;
					float mis_weight = power_heuristic(prev_pdf, pdf_light);

					//float inv_pdf_light = (light_area_pdf.x * cos_light) * light_area_pdf.y / length2;
					//float mis_weight =  power_heuristic_inv_parameter()

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

				//vec3 light_normal = texelFetch(lights_tex, 6 * light_index + 3).xyz;
				
				Light4 l = lights[light_index];

				vec3 light_normal = l.n.xyz;
				//vec3 light_normal = lights[6 * light_index + 3].xyz;

				//neu de cos_light nay se gay ra nhieu vung firefly
				//float cos_light = abs((dot(light_direction, light_normal)));

				//float cos_light = (-1.0f * (dot(light_direction, light_normal)));

				float cos_light = ((dot(light_direction, light_normal)));

                if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_shadow_dfs_final(light_ray, length - eps))
				//if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_shadow_dfs(light_ray, length - eps))
				//if(cos_mtl > 0.0f && cos_light < 0.0f && !hit_shadow(light_ray, length - eps))
				//if(!hit_shadow(light_ray, length - eps))
				{	
					Light4 l = lights[light_index];
					//vec3 light_emission = texelFetch(lights_tex, 6 * light_index + 4).xyz;
					
					vec3 light_emission = l.e.xyz;
					//vec3 light_emission = lights[6 * light_index + 4].xyz;
					
					vec3 light_area_pdf = l.area_pdf.xyz;
					//vec3 light_area_pdf = texelFetch(lights_tex, 6 * light_index + 5).xyz;
					//vec3 light_area_pdf = lights[6 * light_index + 5].xyz;

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

			r.o = new_hit_point;
			r.d = sample_direction;

			//r = Ray(new_hit_point, sample_direction);
		}
		else
			return L;
	}
	return L;
}

void main_halton()
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
	d.x *= tan_fov;//screenResolution.x / screenResolution.y * tan_fov;
	d.y *= tan_fov;

	vec3 rayDir = normalize(d.x * camera.right + d.y * camera.up + camera.forward);

	Ray ray = Ray(camera.position, rayDir);

	//vec3 pixelColor = Path_Tracer2(ray);
	//vec3 pixelColor = Path_Tracer(ray);
	//vec3 pixelColor = brute_force(ray);//Path_Tracer(ray);

	//vec3 pixelColor = path_trace_albedo(ray);

	vec3 accumulate_color = texture(accumulate_tex, tex).xyz;

	vec3 pixelColor = path_trace(ray);

	//vec3 pixelColor = vec3(0,1,1);

	color = pixelColor + accumulate_color;
}

float randf() {
    return fract(sin(dot(seed, vec2(12.9898, 78.233))) * 43758.5453);
}

void main_jitter()
{
	seed = gl_FragCoord.xy + vec2(float(currentFrameSample) * 1.618);

    const float GRID_SIZE = 4.0; 
    
    int sampleIndex = currentFrameSample;
    
    int cellX = sampleIndex % int(GRID_SIZE);
    int cellY = (sampleIndex / int(GRID_SIZE)) % int(GRID_SIZE);
    
    float r1 = rand();
    float r2 = rand();
    
    vec2 subPixelJitter;
    subPixelJitter.x = (float(cellX) + r1) / GRID_SIZE; // Giá trị trả về chuẩn xác trong [0, 1)
    subPixelJitter.y = (float(cellY) + r2) / GRID_SIZE;

   
    vec2 jitter;
    float tx = subPixelJitter.x * 2.0;
    float ty = subPixelJitter.y * 2.0;
    jitter.x = tx < 1.0 ? sqrt(tx) - 1.0 : 1.0 - sqrt(2.0 - tx);
    jitter.y = ty < 1.0 ? sqrt(ty) - 1.0 : 1.0 - sqrt(2.0 - ty);

    
    float aspect = screenResolution.x / screenResolution.y;
    jitter.x *= aspect;

    
    jitter /= (screenResolution * 0.5);

	vec2 d = (2.0f * tex - 1.0f) + jitter;

	float tan_fov = tan(camera.fov * 0.5f);
	d.x *= tan_fov;//screenResolution.x / screenResolution.y * tan_fov;
	d.y *= tan_fov;

	vec3 rayDir = normalize(d.x * camera.right + d.y * camera.up + camera.forward);

	Ray ray = Ray(camera.position, rayDir);

	//vec3 pixelColor = Path_Tracer2(ray);
	//vec3 pixelColor = Path_Tracer(ray);
	//vec3 pixelColor = brute_force(ray);//Path_Tracer(ray);

	//vec3 pixelColor = path_trace_albedo(ray);

	vec3 accumulate_color = texture(accumulate_tex, tex).xyz;

	vec3 pixelColor = path_trace(ray);

	//vec3 pixelColor = vec3(0,1,1);

	color = pixelColor + accumulate_color;
}

float interleavedGradientNoise(vec2 uv, int frameIndex) {
    // Công thức ma trận đặc biệt triệt tiêu tần số thấp, biến khối vuông thành hạt nhiễu động dạng cát rất mịn
    float f = frameIndex;
    uv += vec2(f * 5.588238, f * 3.32512);
    return fract(52.9829189 * fract(dot(uv, vec2(0.0605, 0.0598))));
}

void main_interleave()
{
	int frameIndex = currentFrameSample;
    
    // Lấy mẫu r1, r2 từ tọa độ màn hình thô kết hợp biến tăng Frame
    float r1 = interleavedGradientNoise(gl_FragCoord.xy, frameIndex) * 2.0;
    float r2 = interleavedGradientNoise(gl_FragCoord.xy + vec2(47.12, 93.45), frameIndex) * 2.0;

    vec2 jitter;
    jitter.x = r1 < 1.0 ? sqrt(r1) - 1.0 : 1.0 - sqrt(2.0 - r1);
    jitter.y = r2 < 1.0 ? sqrt(r2) - 1.0 : 1.0 - sqrt(2.0 - r2);

    float aspect = screenResolution.x / screenResolution.y;
    jitter.x *= aspect;

    jitter /= (screenResolution * 0.5);

	vec2 d = (2.0f * tex - 1.0f) + jitter;

	float tan_fov = tan(camera.fov * 0.5f);
	d.x *= tan_fov;//screenResolution.x / screenResolution.y * tan_fov;
	d.y *= tan_fov;

	vec3 rayDir = normalize(d.x * camera.right + d.y * camera.up + camera.forward);

	Ray ray = Ray(camera.position, rayDir);

	//vec3 pixelColor = Path_Tracer2(ray);
	//vec3 pixelColor = Path_Tracer(ray);
	//vec3 pixelColor = brute_force(ray);//Path_Tracer(ray);

	//vec3 pixelColor = path_trace_albedo(ray);

	vec3 accumulate_color = texture(accumulate_tex, tex).xyz;

	vec3 pixelColor = path_trace(ray);

	//vec3 pixelColor = vec3(0,1,1);

	color = pixelColor + accumulate_color;
}

void main()
{
	//main_jitter();
	//main_interleave();
	
}