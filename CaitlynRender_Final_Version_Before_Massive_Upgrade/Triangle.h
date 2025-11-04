#ifndef _TRIANGLE_H_
#define _TRIANGLE_H_

#include <glm\glm.hpp>

using namespace glm;

struct Triangle
{
	/*ivec4 v; //xyz: vertex index, w material index
	ivec4 vt;//xyz, w: reserve, may be for alpha opaque
	ivec4 vn;//xyz normal index (w = 1) or normal vector itself (w = 0) 
			 //w = 0: no need to compute just use vn.xyz directly
			 //w = 1: need to compute because vn is just index

	Triangle() {}
	Triangle(ivec4 v_, ivec4 vt_, ivec4 vn_) : v(v_), vt(vt_), vn(vn_) {}*/

	ivec4 v; //xyz: vertex index, w material index
	ivec4 vn;//xyz normal index (w = 1) or normal vector itself (w = 0) 
			 //w = 0: no need to compute just use vn.xyz directly
			 //w = 1: need to compute because vn is just index
	ivec4 vt;//xyz, w: reserve, may be for alpha opaque
	

	Triangle() {}
	Triangle(ivec4 v_, ivec4 vn_, ivec4 vt_) : v(v_), vn(vn_), vt(vt_) {}
};

struct Triangle2
{
	//ivec4 v; //xyz: vertex index, w material index
			
	ivec4 vn;//xyz normal index (w = 1) or normal vector itself (w = 0) 
			 //w = 0: no need to compute just use vn.xyz directly
			 //w = 1: need to compute because vn is just index
	
	ivec4 vt;//xyz, w: material index

	Triangle2() {}
	Triangle2(ivec4 vn_, ivec4 vt_) : vn(vn_), vt(vt_) {}
};

#endif // !_TRIANGLE_H_
