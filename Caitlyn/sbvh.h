#ifndef _SBVH2_H_
#define _SBVH2_H_
#include "Triangle.h"
#include "BBox.h"
//#include "Sort.h"
#include "FlatNode.h"
#include <vector>
#include <queue>
#include <algorithm>

using namespace std;

#define inf 1e20
#define ci 1
#define ct 4

const int num_bin = 256;
const float inv_num_bin = 1.0f / (float)(num_bin);


struct ObjectSplit
{
	float sah = inf;
	int dim;
	int num_left = 0;
	BBox left_bound;
	BBox right_bound;
};

struct SpatialSplit
{
	float sah = inf;
	int dim;
	float pos;
};

struct SpatialBin
{
	BBox box;
	int enter = 0;
	int exit = 0;
};

int clamp(int v, int lo, int hi)
{
	return v < lo ? lo : v > hi ? hi : v;
}
float clamp(float v, float lo, float hi)
{
	return v < lo ? lo : v > hi ? hi : v;
}

ivec3 clamp3i(ivec3 v, ivec3 lo, ivec3 hi)
{
	return ivec3(clamp(v.x, lo.x, hi.x), clamp(v.y, lo.y, hi.y), clamp(v.z, lo.z, hi.z));
}

float min3(float& a, float& b, float& c)
{
	return minf(a, minf(b, c));
}

vec3 lerp(vec3& a, vec3& b, float t)
{
	//return a * (1.0f - t) + b * t;

	return a + t * (b - a);
}

struct parent_child_node
{
	node* parent = new node();
	NodeSpec spec;

	parent_child_node() { parent = new node(); }
	parent_child_node(node* p, NodeSpec s) : parent(p), spec(s) {}

	//parent_child_node(node* p, node* c, NodeSpec s) : parent(p), child(c), spec(s) {}
};

struct SBVH
{
	vector<Reference> references;

	vector<int> triangle_indices;

	vector<BBox> m_rightBound;

	vector<FlatNode> flat_nodes;

	SpatialBin m_bins[3][num_bin];

	int leaf = 2;

	float m_min_overlap;
	float split_alpha = 0.00001;//1e-5;

	SBVH() {}
	SBVH(vector<Triangle>& trs, vector<vec3>& vertices)
	{
		NodeSpec root_spec;

		int s = trs.size();

		root_spec.num_ref = s;

		references.resize(s);

		for (int i = 0; i < s; ++i)
		{
			references[i].id = i;
			ivec4 v = trs[i].v;
			
			for (int j = 0; j < 3; ++j)
				references[i].box.expand(vertices[v[j]]);
			
			root_spec.box.expand(references[i].box);
		}

		m_min_overlap = root_spec.box.area() * split_alpha;

		//leaf = s > 1000 ? leaf : leaf + 1;

		m_rightBound.resize(maxf(root_spec.num_ref, num_bin) - 1);

		node* root;

		build_iterative(root, root_spec, trs, vertices);

		vector<Triangle> new_trs;

		int new_s = triangle_indices.size();

		new_trs.resize(new_s);

		for (int i = 0; i < new_s; ++i)
			new_trs[i] = trs[triangle_indices[i]];

		trs = new_trs;

		vector<Triangle>().swap(new_trs);

		convert_to_bvh1(root, trs, vertices);

		//int c = 0;
		//count_number_of_leaf_greater_than_2(root, c);

		//cout << "Node greater than 1: " << c << "\n";

		flat_bvh(root);

		//delete_bvh(root);
	}

	/*void count_number_of_leaf_greater_than_2(node* n, int& count)
	{
		if (!n)
			return;
		if (n->leaf)
		{
			if (n->range > 1)
				++count;
			return;
		}

		count_number_of_leaf_greater_than_2(n->nodes[0], count);
		count_number_of_leaf_greater_than_2(n->nodes[1], count);
	}*/

	node* create_leaf(NodeSpec& spec)
	{
		node* n = new node();

		vector<int>& tris = triangle_indices;

		for (int i = 0; i < spec.num_ref; ++i)
		{
			tris.emplace_back(references.back().id);
			references.pop_back();
		}

		n->box = spec.box;
		n->start = tris.size() - spec.num_ref;
		n->range = spec.num_ref;
		n->leaf = 1;

		return n;
	}

	void create_leaf_iterative(node*& n, NodeSpec& spec)
	{
		vector<int>& tris = triangle_indices;

		for (int i = 0; i < spec.num_ref; ++i)
		{
			tris.emplace_back(references.back().id);
			references.pop_back();
		}

		n->box = spec.box;
		n->start = tris.size() - spec.num_ref;
		n->range = spec.num_ref;
		n->leaf = 1;

	}

	int count_leaf()
	{
		int count = 0;
		for (auto& v : flat_nodes)
			if (v.is_leaf())
				++count;
		cout << "Nodes Size: " << flat_nodes.size() << "\n";
		cout << "Leaf Count: " << count << "\n";
		return count;
	}

	void build_iterative(node*& n, NodeSpec& sp, vector<Triangle>& trs, vector<vec3>& vertices)
	{
		n = new node();
		n->box = sp.box;
		n->leaf = 0;

		parent_child_node stk[64];
		int ind = 0;

		//node* dummy = new node();
		stk[ind] = parent_child_node(n, sp);

		while (ind >= 0)
		{
			//cout << ind << "\n";
			auto top = stk[ind--];

			node*& parent = top.parent;


			NodeSpec& top_spec = top.spec;
			
			parent->box = top_spec.box;
			parent->leaf = 0;

			if (top_spec.num_ref <= leaf)
			{
				
				create_leaf_iterative(parent, top_spec);
				continue;
			}

			float node_area = top_spec.box.area();
			float leaf_sah = node_area * top_spec.num_ref;
			float node_sah = 2.0f * node_area;


			ObjectSplit object = find_object_split(top_spec, node_sah);


			BBox overlap = object.left_bound;
			overlap.intersect(object.right_bound);

			SpatialSplit spatial;
			if (overlap.area() >= m_min_overlap)
				spatial = find_spatial_split(top_spec, node_sah, trs, vertices);

			float min_sah = minf(object.sah, spatial.sah);

			NodeSpec left, right;
			if (min_sah == spatial.sah)
				perform_spatial_split(left, right, top_spec, spatial, trs, vertices);			
			else
				perform_object_split(left, right, top_spec, object);

			

			parent->nodes[0] = new node();
			parent->nodes[1] = new node();
			
			stk[++ind] = { parent->nodes[0], left };
			stk[++ind] = { parent->nodes[1], right };
			
		}
		
	}

	void convert_to_bvh1(node*& n, vector<Triangle>& new_trs, vector<vec3>& vertices)
	{
		if (!n)
			return;
		if (n->leaf && n->range == 2)
		{
			n->leaf = 0;

			int left_triangle = n->start;
			int right_triangle = n->start + 1;
			
			ivec4 left_ind = new_trs[left_triangle].v;
			ivec4 right_ind = new_trs[right_triangle].v;

			BBox left_box;
			BBox right_box;

			for (int i = 0; i < 3; ++i)
			{
				left_box.expand(vertices[left_ind[i]]);
				right_box.expand(vertices[right_ind[i]]);
			}

			node* left = new node(left_triangle, 1);
			node* right = new node(right_triangle, 1);

			left->box = left_box;
			right->box = right_box;

			left->leaf = 1;
			right->leaf = 1;

			n->nodes[0] = left;
			n->nodes[1] = right;
			return;
		}

		convert_to_bvh1(n->nodes[0], new_trs, vertices);
		convert_to_bvh1(n->nodes[1], new_trs, vertices);
	}

	void sort_reference(int& num_ref, int& dim)
	{
		sort(references.begin() + references.size() - num_ref, references.end(),
			[dim](const Reference& r1, const Reference& r2)
		{
			float ca = r1.box.c()[dim];
			float cb = r2.box.c()[dim];

			return (ca < cb) || (ca == cb && r1.id < r2.id);
		});
	}

	ObjectSplit find_object_split(NodeSpec& spec, float node_sah)
	{
		ObjectSplit split;

		//Reference* ref = references.data() + references.size() - spec.num_ref;

		//float inv_a = 1.0f / (spec.box.area());
		int start = references.size() - spec.num_ref;
		for (int dim = 0; dim < 3; ++dim)
		{
			sort_reference(spec.num_ref, dim);

			//quick_sort_iterative(references, references.size() - spec.num_ref, references.size(), dim);

			BBox right_box;

			for (int i = spec.num_ref - 1; i > 0; --i)
			{
				right_box.expand(references[start + i].box);
				m_rightBound[i - 1] = right_box;
			}

			BBox left_box;

			for (int i = 1; i < spec.num_ref; i++)
			{
				left_box.expand(references[start + i - 1].box);

				float sah = node_sah + (left_box.area() * i + m_rightBound[i - 1].area() * (spec.num_ref - i));// *inv_a;
				if (sah < split.sah)
				{
					split.sah = sah;
					split.dim = dim;
					split.num_left = i;
					split.left_bound = left_box;
					split.right_bound = m_rightBound[i - 1];
				}
			}
		}
		return split;
	}
	void perform_object_split(NodeSpec& left, NodeSpec& right, NodeSpec& spec, ObjectSplit& split)
	{
		int dim = split.dim;
		
		sort_reference(spec.num_ref, dim);
		
		left.num_ref = split.num_left;
		left.box = split.left_bound;
		right.num_ref = spec.num_ref - split.num_left;
		right.box = split.right_bound;
	}

	void split_reference(Reference& left, Reference& right, Reference& ref, int dim, float pos, vector<Triangle>& trs, vector<vec3>& vertices)
	{
		left.id = right.id = ref.id;
		left.box = right.box = BBox();

		ivec4 v = trs[ref.id].v;

		for (int i = 0; i < 3; ++i)
		{
			vec3 v0 = vertices[v[i]];
			vec3 v1 = vertices[v[(i + 1) % 3]];

			float v0p = v0[dim];
			float v1p = v1[dim];

			if (v0p <= pos)
				left.box.expand(v0);
			if (v0p >= pos)
				right.box.expand(v0);

			if ((v0p < pos && v1p > pos) || (v0p > pos && v1p < pos))
			{
				vec3 t = lerp(v0, v1, clamp((pos - v0p) / (v1p - v0p), 0.0f, 1.0f));
				left.box.expand(t);
				right.box.expand(t);
			}
		}
		left.box.bbox[1][dim] = pos;
		right.box.bbox[0][dim] = pos;
		left.box.intersect(ref.box);
		right.box.intersect(ref.box);
	}

	SpatialSplit find_spatial_split(NodeSpec& spec, float& node_sah, vector<Triangle>& trs, vector<vec3>& vertices)
	{
		vec3 origin = spec.box.bbox[0];
		vec3 bin_size = (spec.box.bbox[1] - origin) * (1.0f / (float)num_bin);//inv_num_bin;
		vec3 inv_bin_size = vec3(1.0f / bin_size.x, 1.0f / bin_size.y, 1.0f / bin_size.z);

		for (int dim = 0; dim < 3; ++dim)
		{
			for (int i = 0; i < num_bin; ++i)
			{
				SpatialBin& bin = m_bins[dim][i];
				bin.box = BBox();
				bin.enter = 0;
				bin.exit = 0;
			}
		}

		for (int refIdx = references.size() - spec.num_ref; refIdx < references.size(); ++refIdx)
		{
			const Reference& ref = references[refIdx];

			ivec3 first_bin = clamp3i(ivec3((ref.box.bbox[0] - origin) * inv_bin_size), ivec3(0), ivec3(num_bin - 1));
			ivec3 last_bin = clamp3i(ivec3((ref.box.bbox[1] - origin) * inv_bin_size), first_bin, ivec3(num_bin - 1));

			for (int dim = 0; dim < 3; ++dim)
			{
				Reference curr_ref = ref;
				for (int i = first_bin[dim]; i < last_bin[dim]; ++i)
				{
					Reference left_ref, right_ref;
					split_reference(left_ref, right_ref, curr_ref, dim, origin[dim] + bin_size[dim] * (float)(i + 1), trs, vertices);
					m_bins[dim][i].box.expand(left_ref.box);
					curr_ref = right_ref;
				}
				m_bins[dim][last_bin[dim]].box.expand(curr_ref.box);
				m_bins[dim][first_bin[dim]].enter++;
				m_bins[dim][last_bin[dim]].exit++;
			}
		}
		SpatialSplit split;

		for (int dim = 0; dim < 3; ++dim)
		{
			BBox right_box;
			for (int i = num_bin - 1; i > 0; --i)
			{
				right_box.expand(m_bins[dim][i].box);
				m_rightBound[i - 1] = right_box;
			}

			BBox left_box;
			int left_num = 0;
			int right_num = spec.num_ref;

			for (int i = 1; i < num_bin; ++i)
			{
				left_box.expand(m_bins[dim][i - 1].box);
				left_num += m_bins[dim][i - 1].enter;
				right_num -= m_bins[dim][i - 1].exit;

				float sah = node_sah + left_box.area() * left_num + m_rightBound[i - 1].area() * right_num;

				if (sah < split.sah)
				{
					split.sah = sah;
					split.dim = dim;
					split.pos = origin[dim] + bin_size[dim] * (float)i;
				}
			}
		}
		return split;
	}

	void perform_spatial_split(NodeSpec& left, NodeSpec& right, NodeSpec& spec, SpatialSplit& split, vector<Triangle>& trs, vector<vec3>& vertices)
	{
		vector<Reference>& refs = references;

		int left_start = refs.size() - spec.num_ref;
		int left_end = left_start;
		int right_start = refs.size();
		left.box = right.box = BBox();

		int dim = split.dim;
		float pos = split.pos;

		for (int i = left_end; i < right_start; ++i)
		{
			if (refs[i].box.bbox[1][dim] <= pos)
			{
				left.box.expand(refs[i].box);
				swap(refs[i], refs[left_end++]);
			}
			else if (refs[i].box.bbox[0][dim] >= pos)
			{
				right.box.expand(refs[i].box);
				swap(refs[i--], refs[--right_start]);
			}
		}

		while (left_end < right_start)
		{
			Reference lref, rref;
			split_reference(lref, rref, refs[left_end], dim, pos, trs, vertices);

			BBox lub = left.box;
			BBox rub = right.box;
			BBox ldb = left.box;
			BBox rdb = right.box;

			lub.expand(refs[left_end].box);
			rub.expand(refs[left_end].box);
			ldb.expand(lref.box);
			rdb.expand(rref.box);

			float lac = left_end - left_start;
			float rac = refs.size() - right_start;
			float lbc = lac + 1;
			float rbc = rac + 1;

			float unsplit_left_sah = lub.area() * lbc + right.box.area() * rac;
			float unsplit_right_sah = left.box.area() * lac + rub.area() * rbc;
			float duplicate_sah = ldb.area() * lbc + rdb.area() * rbc;

			float min_sah = min3(unsplit_left_sah, unsplit_right_sah, duplicate_sah);

			if (min_sah == unsplit_left_sah)
			{
				left.box = lub;
				left_end++;
			}
			else if (min_sah == unsplit_right_sah)
			{
				right.box = rub;
				swap(refs[left_end], refs[--right_start]);
			}
			else
			{
				left.box = ldb;
				right.box = rdb;
				refs[left_end++] = lref;
				refs.emplace_back(rref);
			}
		}
		left.num_ref = left_end - left_start;
		right.num_ref = refs.size() - right_start;
	}
	void flat_bvh(node*& n)
	{
		queue<node*> queue_node;

		queue_node.emplace(n);

		int left_index = 0;

		while (!queue_node.empty())
		{
			node* front = queue_node.front();

			queue_node.pop();

			if (front->leaf)
			{
				FlatNode flatnode;

				flatnode.box_min = vec4(front->box.bbox[0], front->start);
				flatnode.box_max = vec4(front->box.bbox[1], front->range);

				flat_nodes.emplace_back(flatnode);
			}
			else
			{
				FlatNode flatnode;

				flatnode.box_min = vec4(front->box.bbox[0], ++left_index);
				flatnode.box_max = vec4(front->box.bbox[1], 0);

				flat_nodes.emplace_back(flatnode);

				queue_node.emplace(front->nodes[0]);
				queue_node.emplace(front->nodes[1]);

				++left_index;
			}
		}
		//cout << "Completely FLat !\n";
	}

	void clean_up()
	{
		vector<Reference>().swap(references);

		vector<int>().swap(triangle_indices);

		vector<BBox>().swap(m_rightBound);

		vector<FlatNode>().swap(flat_nodes);
	}
	void delete_bvh(node*& n)
	{
		if (n == NULL)
			return;
		delete_bvh(n->nodes[0]);
		delete_bvh(n->nodes[1]);

		delete(n);
	}
};

#endif // !_SBVH2_H_

