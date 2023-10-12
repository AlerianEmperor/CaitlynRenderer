#ifndef _FLAT_NODE_H_
#define _FLAT_NODE_H_
#include <glm\glm.hpp>

using namespace glm;

struct Reference
{
	int id;
	BBox box;
};

struct NodeSpec
{
	int num_ref = 0;
	BBox box;
};

struct node
{
	BBox box;

	node *nodes[2] = { NULL, NULL };
	unsigned start : 21;
	unsigned range : 8;
	unsigned axis : 2;
	unsigned leaf : 1;


	node() { start = 0, range = 0; }
	node(int s, int r) : start(s), range(r) {}
};

struct FlatNode
{
	vec4 box_min;//start
	vec4 box_max;//range

	FlatNode() {}
	FlatNode(vec4 bmin, vec4 bmax) : box_min(bmin), box_max(bmax) {}
	float area()
	{
		//no need to multiply by 2, because if everything is not multiply by 2, the result is the same
		vec4 d = box_max - box_min;

		return d.x * (d.y + d.z) + d.y * d.z;
	}
	vec3 bmin()
	{
		return vec3(box_min.x, box_min.y, box_min.z);
	}
	vec3 bmax()
	{
		return vec3(box_max.x, box_max.y, box_max.z);
	}
	vec3 center()
	{
		return vec3(box_min.x + box_max.x, box_min.y + box_max.y, box_min.z + box_max.z) * 0.5f;
	}
	bool is_leaf()
	{
		return box_max.w != 0;
	}
	BBox box()
	{
		vec3 bmin(box_min.x, box_min.y, box_min.z);
		vec3 bmax(box_max.x, box_max.y, box_max.z);

		return BBox(bmin, bmax);
	}
};

FlatNode unify_box(FlatNode& n1, FlatNode& n2)
{
	vec3 bmin = vec3(maxf(n1.box_min.x, n2.box_min.x), maxf(n1.box_min.y, n2.box_min.y), maxf(n1.box_min.z, n2.box_min.z));
	vec3 bmax = vec3(minf(n1.box_max.x, n2.box_max.x), maxf(n1.box_max.y, n2.box_max.y), maxf(n1.box_max.z, n2.box_max.z));

	return FlatNode(vec4(bmin, 0), vec4(bmax, 0));
}

#endif // !_FLAT_NODE_H

