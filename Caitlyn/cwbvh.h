#ifndef _CWBVH_H_
#define _CWBVH_H_
#include "sbvh.h"
#include <iostream>

using namespace std;

typedef unsigned char byte;
char INVALID = '-';

struct node8
{
	vec3 p;
	byte e[3];
	byte imask;

	unsigned int child_base_index;
	unsigned int triangle_base_index;

	byte meta[8] = {};

	byte quantized_min_x[8] = {}, quantized_max_x[8] = {};
	byte quantized_min_y[8] = {}, quantized_max_y[8] = {};
	byte quantized_min_z[8] = {}, quantized_max_z[8] = {};

	bool is_leaf(int child_ind)
	{
		//3 highgest bit may be : 001, 011, 111, equal to 1, 2 and 3 triangle
		//only inner node have bit 3 and 4 set
		//which leaves the last three low bit for leaf node, range from : 001 to 111

		//so (meta[child_ind] & 0b00011111) < 8; is enough
		//24 is overkill since only inner node have 3rd and 4th bit set: 11000

		return (meta[child_ind] & 0b00011111) < 24;
	}
};

enum node_type : char { LEAF, INTERNAL, DISTRIBUTE };

struct Decision
{
	node_type type;

	char distribute_left;
	char distribute_right;

	float cost;
};

struct CWBVH
{
	vector<node8> nodes;
	vector<int> triangle_indices;
	vector<Decision> decisions;

	CWBVH() {}
	void convert(SBVH& bvh)
	{
		vector<node8>().swap(nodes);
		nodes.reserve(bvh.flat_nodes.size());
		
		vector<int>().swap(triangle_indices);
		triangle_indices.reserve(bvh.triangle_indices.size());

		decisions.resize(bvh.flat_nodes.size() * 7);

		calculate_cost(0, bvh);

		collapse(bvh, 0, 0);
		
		nodes.emplace_back();
	}
	//checked
	int calculate_cost(int node_ind, SBVH& bvh)
	{
		FlatNode& node = bvh.flat_nodes[node_ind];
		int num_primitive;

		if (node.is_leaf())
		{
			num_primitive = node.box_max.w;//equal 1 any way!

			float cost_leaf = node.area() * float(num_primitive);

			for (int i = 0; i < 7; ++i)
			{
				decisions[7 * node_ind + i].type = LEAF;
				decisions[7 * node_ind + i].cost = cost_leaf;
			}
		}
		else
		{
			int left = node.box_min.w;
			int right = left + 1;

			num_primitive = calculate_cost(left, bvh) + calculate_cost(right, bvh);

			//i = 0

			float cost_leaf = num_primitive <= 3 ? float(num_primitive) * node.area() : inf;

			float cost_distribute = inf;

			char distribute_left = -1;//INVALID;
			char distribute_right = -1;// INVALID;

			for (int k = 0; k < 7; ++k)
			{
				float c = decisions[7 * left + k].cost + decisions[7 * right + 6 - k].cost;

				if (c < cost_distribute)
				{
					cost_distribute = c;

					distribute_left = k;
					distribute_right = 6 - k;
				}
			}

			float cost_internal = cost_distribute + node.area();

			if (cost_leaf < cost_internal)
			{
				decisions[7 * node_ind].type = LEAF;
				decisions[7 * node_ind].cost = cost_leaf;
			}
			else
			{
				decisions[7 * node_ind].type = INTERNAL;
				decisions[7 * node_ind].cost = cost_internal;
			}

			decisions[7 * node_ind].distribute_left = distribute_left;
			decisions[7 * node_ind].distribute_right = distribute_right;

			//i = 1..7
			for (int i = 1; i < 7; ++i)
			{
				float cost_distribute = decisions[7 * node_ind + i - 1].cost;

				char distribute_left = -1;//INVALID;
				char distribute_right = -1;//INVALID;

				for (int k = 0; k < i; ++k)
				{
					float c = decisions[7 * left + k].cost +
						      decisions[7 * right + i - k - 1].cost;

					if (c < cost_distribute)
					{
						cost_distribute = c;

						distribute_left = k;
						distribute_right = i - k - 1;
					}
				}

				decisions[7 * node_ind + i].cost = cost_distribute;

				if (distribute_left != INVALID)
				{
					decisions[7 * node_ind + i].type = DISTRIBUTE;
					decisions[7 * node_ind + i].distribute_left = distribute_left;
					decisions[7 * node_ind + i].distribute_right = distribute_right;
				}
				else
					decisions[7 * node_ind + i] = decisions[7 * node_ind + i - 1];
			}
		}
		
		return num_primitive;
	}

	void get_children(int node_ind, SBVH& bvh, int children[8], int& child_count, int i)
	{
		FlatNode& node = bvh.flat_nodes[node_ind];

		if (node.is_leaf())
		{
			children[child_count++] = node_ind;
			return;
		}

		char distribute_left = decisions[7 * node_ind + i].distribute_left;
		char distribute_right = decisions[7 * node_ind + i].distribute_right;

		//if (distribute_left < 0 || distribute_left >= 7 || distribute_right < 0 || distribute_right >= 7)
		//	return;

		int left = node.box_min.w;
		
		if (decisions[7 * left + distribute_left].type == DISTRIBUTE)
			get_children(left, bvh, children, child_count, distribute_left);
		else
			children[child_count++] = left;

		int right = left + 1;

		if (decisions[7 * right + distribute_right].type == DISTRIBUTE)
			get_children(right, bvh, children, child_count, distribute_right);
		else
			children[child_count++] = right;
	}
	//xong
	void order_children(int node_ind, SBVH& bvh, int children[8], int child_count)
	{
		vec3 p = bvh.flat_nodes[node_ind].center();

		float cost[8][8] = {};

		for (int c = 0; c < child_count; c++)
		{
			for (int s = 0; s < 8; s++)
			{
				vec3 direction(
					(s & 0b100) ? -1.0f : 1.0f,
					(s & 0b010) ? -1.0f : 1.0f,
					(s & 0b001) ? -1.0f : 1.0f
				);

				cost[c][s] = dot(bvh.flat_nodes[children[c]].box().c() - p, direction);
			}
		}

		int assignment[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };//{ INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID };
		bool slot_filled[8] = {};

		while (true)
		{
			float min_cost = inf;

			int min_slot = -1;//INVALID;
			int min_index = -1;// INVALID;

			for (int c = 0; c < child_count; ++c)
			{
				if (assignment[c] == -1)//INVALID)
				{
					for (int s = 0; s < 8; ++s)
					{
						if (!slot_filled[s] && cost[c][s] < min_cost)
						{
							min_cost = cost[c][s];

							min_slot = s;
							min_index = c;
						}
					}
				}
			}

			if (min_slot == -1)//INVALID)
				break;

			slot_filled[min_slot] = true;
			assignment[min_index] - min_slot;
		}
		int new_children[8];

		for (int i = 0; i < 8; ++i)
		{
			new_children[i] = children[i];
			children[i] = -1;//INVALID;
		}

		for (int i = 0; i < child_count; ++i)
		{
			//if (assignment[i] != INVALID && new_children[i] != INVALID)
			children[assignment[i]] = new_children[i];
		}
	}
	//xong
	int count_primitive(int node_ind, SBVH& bvh)
	{
		FlatNode& node = bvh.flat_nodes[node_ind];

		if (node.is_leaf())
		{
			int start = node.box_min.w;
			int range = node.box_max.w;

			for (unsigned i = 0; i < range; ++i)
				triangle_indices.emplace_back(bvh.triangle_indices[start + i]);

			return range;
		}

		int left = node.box_min.w;
		return count_primitive(left, bvh) +
			   count_primitive(left + 1, bvh);
	}
	//xong
	void collapse(SBVH& bvh, int node_ind_bvh2, int node_ind_bvh8)
	{
		node8& node = nodes[node_ind_bvh8];

		BBox& box = bvh.flat_nodes[node_ind_bvh2].box();

		node.p = box.bbox[0];

		constexpr int Nq = 8;
		constexpr float denom = 1.0f / float((1 << Nq) - 1);

		float ex = exp2f(ceilf(log2f((box.bbox[1].x - box.bbox[0].x) * denom)));
		float ey = exp2f(ceilf(log2f((box.bbox[1].y - box.bbox[0].y) * denom)));
		float ez = exp2f(ceilf(log2f((box.bbox[1].z - box.bbox[0].z) * denom)));

		vec3 e(ex, ey, ez);

		vec3 inv_e(1.0f / ex, 1.0f / ey, 1.0f / ez);

		unsigned u_ex = (unsigned)(e.x);
		unsigned u_ey = (unsigned)(e.y);
		unsigned u_ez = (unsigned)(e.z);

		//quantization

		node.e[0] = u_ex >> 23;
		node.e[1] = u_ey >> 23;
		node.e[2] = u_ez >> 23;
		
		int child_count = 0;
		int children[8] = { -1,-1,-1,-1,-1,-1,-1,-1 };//{ INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID, INVALID };

		get_children(node_ind_bvh2, bvh, children, child_count, 0);

		order_children(node_ind_bvh2, bvh, children, child_count);

		node.imask = 0;

		node.triangle_base_index = unsigned(triangle_indices.size());
		node.child_base_index = unsigned(nodes.size());

		int num_internal_nodes = 0;
		int num_triangles = 0;

		for (int i = 0; i < 8; ++i)
		{
			int child_index = children[i];

			if (child_index == -1)//INVALID)
				continue;

			BBox& child_bbox = bvh.flat_nodes[child_index].box();

			node.quantized_min_x[i] = byte(floorf((child_bbox.bbox[0].x - node.p.x) * inv_e.x));
			node.quantized_min_y[i] = byte(floorf((child_bbox.bbox[0].y - node.p.y) * inv_e.y));
			node.quantized_min_z[i] = byte(floorf((child_bbox.bbox[0].z - node.p.z) * inv_e.z));

			node.quantized_max_x[i] = byte(floorf((child_bbox.bbox[1].x - node.p.x) * inv_e.x));
			node.quantized_max_y[i] = byte(floorf((child_bbox.bbox[1].y - node.p.y) * inv_e.y));
			node.quantized_max_z[i] = byte(floorf((child_bbox.bbox[1].z - node.p.z) * inv_e.z));

			switch (decisions[7 * child_index].type)
			{
				case LEAF:
				{
					int triangle_count = count_primitive(child_index, bvh);

					//ASSERT triangle_count

					/*if (triangle_count <= 0 || triangle_count > 3)
					{
						cout << "Leaf is too big!\n";
						return;
					}*/

					//3 highest is unary count of number of triangle
					//001: 1 triangle
					//011: 2 triangles
					//111: 3 triangles
					//maximum 3 triangle per leaf is allowed

					for (int j = 0; j < triangle_count; ++j)
						node.meta[i] |= (1 << (j + 5));

					node.meta[i] |= num_triangles;

					num_triangles += triangle_count;
					break;
				}
				case INTERNAL:
				{
					node.meta[i] = (i + 24) | 0b00100000;

					node.imask |= (1 << i);
					num_internal_nodes++;
					break;
				}
			}

			for (int i = 0; i < num_internal_nodes; ++i)
				nodes.emplace_back();

			node = nodes[node_ind_bvh8];

			int offset = 0;

			for (int i = 0; i < 8; ++i)
			{
				int child_index = children[i];

				if (child_index == -1)//INVALID)
					continue;

				if (node.imask & (1 << i))
					collapse(bvh, node.child_base_index + offset++, child_index);
			}
		}
	}
};

#endif // !_CWBVH_H_

