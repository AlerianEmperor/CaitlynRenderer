#ifndef _SBVH_OPTIMIZE_H_
#define _SBVH_OPTIMIZE_H_
#include "Triangle.h"
#include "BBox.h"
//#include "Sort.h"
#include "FlatNode.h"
#include <vector>
#include <queue>
#include <algorithm>
#include <stack>

using namespace std;

#define inf 1e20
#define ci 1
#define ct 4

const int num_bin = 256;//256
const float inv_num_bin = 1.0f / (float)(num_bin);

bool use_soa = false;

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

struct BVH_SOA
{
	vector<float> min_x, min_y, min_z;
	vector<float> max_x, max_y, max_z;

	vector<float> start, range;
	//vector<int> start;
	//vector<int> range;
};

struct SBVH
{
	vector<Reference> references;

	vector<int> triangle_indices;

	vector<BBox> m_rightBound;

	//if use bvh2
	vector<FlatNode> flat_nodes;

	//if use bvh4
	vector<FlatNode4> flat_nodes4;

	SpatialBin m_bins[3][num_bin];

	BVH_SOA bvh_soa;

	int leaf = 2;

	float m_min_overlap;
	float split_alpha = 0.00001;//1e-5;

	int nodeCount;
	int tree_height;

	bool use_bvh2 = true;
	
	SBVH() {}
	SBVH(vector<Triangle>& trs, vector<vec3>& vertices, bool use_bvh2_ = true)//false = wbvh
	{
		use_bvh2 = use_bvh2_;

		NodeSpec root_spec;

		//leaf = std::max(2, (int)std::min(32.0, sqrt((double)trs.size()) * 0.5));

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

		//convert_to_bvh1(root, trs, vertices);

		//flat_bvh(root);
		//int index_counter = 0;
		//flat_bvh_dfs(root, index_counter);

		//flat_bvh_dfs(root);

		//tot
		if (use_bvh2) 
			flat_bvh_dfs4(root);
		else
			//flat_bvh_to_bvh4(root);
			flat_bvh_to_bvh4_dfs(root);

		//flat_bvh_dfs_stackless2(root);

		//flat_bvh_dfs_parent_pointers(root);

		//flat_bvh_dfs5(root);

		//flat_bvh_dfs4_stackless(root);

		//flat_bvh_dfs4_center_extend(root);

		//flat_bvh_adjacent_children(root);

		delete_bvh(root);
	}

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

		tree_height = height(n);

		cout << "tree height " << tree_height << "\n";


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
		nodeCount = flat_nodes.size();
		cout << "nodeCount: " << nodeCount << "\n";

		if (use_soa)
			convert_to_soa();
		//cout << "Completely FLat !\n";
	}

	// Thay th? hàm flat_bvh c? c?a b?n b?ng hàm ?? quy DFS này:
	void flat_bvh_dfs(node* n, int& current_index)
	{
		if (!n) return;

		int my_index = current_index;
		FlatNode flatnode;

		if (n->leaf)
		{
			flatnode.box_min = vec4(n->box.bbox[0], n->start);
			flatnode.box_max = vec4(n->box.bbox[1], n->range);
			flat_nodes.push_back(flatnode);
		}
		else
		{
			// Gán t?m th?i, tí n?a con trái s? b?t ??u ngay ? v? trí ti?p theo (my_index + 1)
			flatnode.box_min = vec4(n->box.bbox[0], current_index + 1);
			flatnode.box_max = vec4(n->box.bbox[1], 0);
			flat_nodes.push_back(flatnode);

			current_index++;
			flat_bvh_dfs(n->nodes[0], current_index); // Con trái
			//current_index++;
			// Sau khi con trái và toàn b? nhánh c?a nó ???c ??y vào m?ng,
			// v? trí hi?n t?i c?a current_index + 1 chính là v? trí c?a con ph?i!
			// Ta c?p nh?t l?i tr??ng start c?a nút hi?n t?i ?? tr? ?úng
			// Tuy nhiên ?? t?i ?u nh?t: N?u ép con trái luôn c?nh con ph?i, b?n không c?n ??i.

			current_index++;
			flat_bvh_dfs(n->nodes[1], current_index); // Con ph?i
			//current_index++;
		}
	}

	// G?i hàm này trong hàm kh?i t?o thay cho flat_bvh(root) c?:
	// int index_counter = 0;
	// flat_bvh_dfs(root, index_counter);

	void flat_bvh_dfs2(node*& n)
	{
		if (!n) return;

		// S? d?ng std::stack thay cho std::queue ?? chuy?n sang thu?t toán DFS
		stack<node*> stack_node;
		stack_node.emplace(n);

		int left_index = 0;

		while (!stack_node.empty())
		{
			node* front = stack_node.top();
			stack_node.pop();

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
				// Trong ki?n trúc DFS ph?ng: Con bên trái s? n?m ngay k? ti?p (??a ch? hi?n t?i + 1)
				// Nên bi?n left_index này s? ??i di?n cho v? trí nh?y c?a con bên PH?I trong m?ng
				left_index += 2;

				flatnode.box_min = vec4(front->box.bbox[0], left_index);
				flatnode.box_max = vec4(front->box.bbox[1], 0);
				flat_nodes.emplace_back(flatnode);

				// ??y con bên Ph?i vào tr??c, con bên Trái vào sau trong stack
				// Quy t?c Stack (LIFO) ??m b?o con bên Trái s? ???c b?c ra x? lý tr??c, 
				// n?m ngay sát c?nh nút cha trong m?ng ph?ng d? li?u
				stack_node.emplace(front->nodes[1]); // Con ph?i
				stack_node.emplace(front->nodes[0]); // Con trái
			}
		}

		nodeCount = flat_nodes.size();
		cout << "nodeCount (DFS Flat): " << nodeCount << "\n";

		if (use_soa)
			convert_to_soa();
	}

	void flat_bvh_dfs3(node* root) {
		if (!root) return;

		struct StackEntry {
			node* ptr;
			int parent_flat_idx;
			int is_right_child;
		};

		std::vector<StackEntry> stk;
		stk.push_back({ root, -1, 0 });

		flat_nodes.clear();

		while (!stk.empty()) {
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int current_idx = flat_nodes.size();

			FlatNode flatnode;
			if (curr->leaf) {
				flatnode.box_min = vec4(curr->box.bbox[0], curr->start);
				// w âm ho?c ?ánh d?u ?? GPU bi?t là node lá, và ch?a s? l??ng tam giác (range)
				flatnode.box_max = vec4(curr->box.bbox[1], curr->range);
			}
			else {
				flatnode.box_min = vec4(curr->box.bbox[0], 0); // S? c?p nh?t index con ph?i sau
				flatnode.box_max = vec4(curr->box.bbox[1], 0);
			}

			flat_nodes.push_back(flatnode);

			// C?p nh?t liên k?t t? nút cha
			if (top.parent_flat_idx != -1) {
				if (top.is_right_child) {
					// N?u ?ây là con ph?i, c?p nh?t thông tin nh?y b??c vào nút cha
					flat_nodes[top.parent_flat_idx].box_max.w = current_idx;
				}
				else {
					// N?u ?ây là con trái, nó luôn b?ng parent_idx + 1 (m?c ??nh theo DFS)
					flat_nodes[top.parent_flat_idx].box_min.w = current_idx;
				}
			}

			if (!curr->leaf) {
				// ??y con ph?i vào tr??c, con trái vào sau ?? con trái ???c x? lý tr??c (DFS)
				stk.push_back({ curr->nodes[1], current_idx, 1 });
				stk.push_back({ curr->nodes[0], current_idx, 0 });
			}
		}
	}

	void flat_bvh_dfs4(node* n)
	{
		if (!n) return;

		flat_nodes.clear();

		// S? d?ng vector làm stack ?? duy?t cây theo chi?u sâu (DFS)
		struct StackEntry {
			node* ptr;
			int parent_flat_idx; // ?? quay l?i c?p nh?t ch? s? con ph?i cho nút cha
		};

		vector<StackEntry> stk;
		stk.push_back({ n, -1 });

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int current_idx = flat_nodes.size(); // V? trí Node ID hi?n t?i trong m?ng

			FlatNode flatnode;
			if (curr->leaf)
			{
				flatnode.box_min = vec4(curr->box.bbox[0], curr->start);
				flatnode.box_max = vec4(curr->box.bbox[1], curr->range); // range > 0 là nút lá
			}
			else
			{
				flatnode.box_min = vec4(curr->box.bbox[0], 0); // S? ghi ?è ch? s? con ph?i vào ?ây sau
				flatnode.box_max = vec4(curr->box.bbox[1], 0); // range = 0 ngh?a là nút trong
			}

			flat_nodes.emplace_back(flatnode);

			// Nút này là con ph?i c?a nút cha tr??c ?ó -> c?p nh?t ch? s? nh?y cho nút cha
			if (top.parent_flat_idx != -1)
			{
				// Vì chúng ta ??y con ph?i vào tr??c, nên khi con ph?i ???c lôi ra kh?i stack,
				// flat_nodes.size() hi?n t?i chính là v? trí b?t ??u c?a nhánh con ph?i!
				flat_nodes[top.parent_flat_idx].box_min.w = current_idx;
			}

			if (!curr->leaf)
			{
				// B?T BU?C: ??y con ph?i vào tr??c, con trái vào sau
				// Nh? v?y khi pop stack, con trái s? ???c x? lý tr??c -> ??t chu?n DFS (Con trái = current_idx + 1)
				stk.push_back({ curr->nodes[1], current_idx }); // Nhánh Ph?i
				stk.push_back({ curr->nodes[0], -1 });          // Nhánh Trái (Không c?n c?p nh?t parent vì nó luôn là idx + 1)
			}
		}

		nodeCount = flat_nodes.size();
		cout << "DFS NodeCount: " << nodeCount << "\n";
	}
	
	void flat_bvh_dfs5(node* root)
	{
		if (!root) return;

		flat_nodes.clear();
		// 1. T?i ?u b? nh?: Tính tr??c t?ng s? node và c?p phát m?t l?n duy nh?t
		int total_nodes = count_nodes(root);
		flat_nodes.resize(total_nodes);

		// 2. Thu?t toán d?ch chuy?n cây không dùng Stack (Thread-safe & Stackless)
		node* curr = root;
		int flat_idx = 0;

		// Con tr? dùng ?? theo dõi nút cha g?n nh?t có con bên ph?i ch?a ???c duy?t
		// Ta t?m th?i m??n m?t tr??ng trong c?u trúc node ho?c qu?n lý b?ng con tr? t?m
		// Tuy nhiên, cách s?ch nh?t là dùng m?ng t?nh (Array-based Stack) v?i kích th??c c? ??nh (64) 
		// thay vì std::vector ?? tri?t tiêu chi phí Heap:
		struct FastStack {
			node* ptr;
			int parent_flat_idx;
		} fast_stk[64];

		int top_ind = 0;
		fast_stk[top_ind] = { root, -1 };

		while (top_ind >= 0)
		{
			auto top = fast_stk[top_ind--];
			curr = top.ptr;
			int current_idx = flat_idx++;

			if (curr->leaf)
			{
				flat_nodes[current_idx].box_min = vec4(curr->box.bbox[0], curr->start);
				flat_nodes[current_idx].box_max = vec4(curr->box.bbox[1], curr->range);
			}
			else
			{
				flat_nodes[current_idx].box_min = vec4(curr->box.bbox[0], 0);
				flat_nodes[current_idx].box_max = vec4(curr->box.bbox[1], 0);
			}

			// C?p nh?t ch? s? nh?y c?a con ph?i cho nút cha (Gi? nguyên logic hoàn h?o c?a b?n)
			if (top.parent_flat_idx != -1)
			{
				flat_nodes[top.parent_flat_idx].box_min.w = current_idx;
			}

			if (!curr->leaf)
			{
				// T?i ?u hóa Cache Hit (Tùy ch?n): Có th? ki?m tra tâm BBox ? ?ây ?? quy?t ??nh
				// xem nên ??y nodes[1] hay nodes[0] vào tr??c nh?m t?i ?u h??ng duy?t tia.

				fast_stk[++top_ind] = { curr->nodes[1], current_idx }; // Nhánh Ph?i
				fast_stk[++top_ind] = { curr->nodes[0], -1 };          // Nhánh Trái
			}
		}

		nodeCount = flat_nodes.size();
		cout << "DFS T?i ?u B? Nh? - NodeCount: " << nodeCount << "\n";

		if (use_soa)
			convert_to_soa();
	}

	void flat_bvh_dfs4_stackless(node* root)
	{
		if (!root) return;

		flat_nodes.clear();
		// C?p phát c?ng dung l??ng tr??c ?? tránh vector b? reallocate liên t?c làm gi?m hi?u n?ng
		flat_nodes.resize(count_nodes(root));

		struct StackEntry {
			node* ptr;
			int skip_idx; // Ch? s? s? nh?y t?i n?u tia TR??T (miss) bounding box c?a nút này
		};

		vector<StackEntry> stk;
		// Nút g?c (Root) n?u tr??t thì thoát vòng l?p, nên skip_idx = -1
		stk.push_back({ root, -1 });

		int flat_idx_counter = 0;

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int my_idx = flat_idx_counter++;

			FlatNode flatnode;
			if (curr->leaf)
			{
				// NÚT LÁ: 
				// box_min.w l?u v? trí start c?a tam giác
				// box_max.w l?u range (>0) ?? GPU bi?t ?ây là nút lá và l?p v? tam giác
				flatnode.box_min = vec4(curr->box.bbox[0], (float)curr->start);
				flatnode.box_max = vec4(curr->box.bbox[1], (float)curr->range);
			}
			else
			{
				// NÚT TRONG (Internal Node):
				// box_min.w l?u skip_index: N?u tr??t nút này, tia nh?y th?ng t?i ?ây (b? qua toàn b? con trái + ph?i)
				// box_max.w = 0.0f ?? ?ánh d?u ?ây là nút trong (không ph?i nút lá)
				flatnode.box_min = vec4(curr->box.bbox[0], (float)top.skip_idx);
				flatnode.box_max = vec4(curr->box.bbox[1], 0.0f);
			}

			flat_nodes[my_idx] = flatnode;

			if (!curr->leaf)
			{
				// T?i ?u hóa tính toán ch? s? nh?y không c?n quay l?i s?a (Backpatching):
				// Trong m?ng ph?ng DFS, con bên ph?i s? luôn b?t ??u n?m ngay sau khi toàn b? nhánh con bên trái k?t thúc.
				int left_subtree_size = count_nodes(curr->nodes[0]);
				int right_child_flat_idx = my_idx + 1 + left_subtree_size;

				// 1. ??y CON PH?I vào tr??c: N?u con ph?i tr??t, m?c tiêu c?a nó là skip_idx c?a cha
				stk.push_back({ curr->nodes[1], top.skip_idx });

				// 2. ??y CON TRÁI vào sau (?? x? lý tr??c): N?u con trái tr??t, m?c tiêu c?a nó là ??u nhánh CON PH?I!
				stk.push_back({ curr->nodes[0], right_child_flat_idx });
			}
		}

		nodeCount = flat_nodes.size();
		cout << "Optimized DFS NodeCount: " << nodeCount << "\n";
	}

	void flat_bvh_dfs4_center_extend(node* n)
	{
		if (!n) return;

		flat_nodes.clear();

		// S? d?ng vector làm stack ?? duy?t cây theo chi?u sâu (DFS)
		struct StackEntry {
			node* ptr;
			int parent_flat_idx; // ?? quay l?i c?p nh?t ch? s? con ph?i cho nút cha
		};

		vector<StackEntry> stk;
		stk.push_back({ n, -1 });

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int current_idx = flat_nodes.size(); // V? trí Node ID hi?n t?i trong m?ng

			FlatNode flatnode;
			if (curr->leaf)
			{
				flatnode.box_min = vec4((curr->box.bbox[1] + curr->box.bbox[0]) * 0.5f, curr->start);
				flatnode.box_max = vec4(curr->box.bbox[1] - curr->box.bbox[0], curr->range); // range > 0 là nút lá
			}
			else
			{
				flatnode.box_min = vec4((curr->box.bbox[1] + curr->box.bbox[0]) * 0.5f, 0); // S? ghi ?è ch? s? con ph?i vào ?ây sau
				flatnode.box_max = vec4(curr->box.bbox[1] - curr->box.bbox[0], 0); // range = 0 ngh?a là nút trong
			}

			flat_nodes.emplace_back(flatnode);

			// Nút này là con ph?i c?a nút cha tr??c ?ó -> c?p nh?t ch? s? nh?y cho nút cha
			if (top.parent_flat_idx != -1)
			{
				// Vì chúng ta ??y con ph?i vào tr??c, nên khi con ph?i ???c lôi ra kh?i stack,
				// flat_nodes.size() hi?n t?i chính là v? trí b?t ??u c?a nhánh con ph?i!
				flat_nodes[top.parent_flat_idx].box_min.w = current_idx;
			}

			if (!curr->leaf)
			{
				// B?T BU?C: ??y con ph?i vào tr??c, con trái vào sau
				// Nh? v?y khi pop stack, con trái s? ???c x? lý tr??c -> ??t chu?n DFS (Con trái = current_idx + 1)
				stk.push_back({ curr->nodes[1], current_idx }); // Nhánh Ph?i
				stk.push_back({ curr->nodes[0], -1 });          // Nhánh Trái (Không c?n c?p nh?t parent vì nó luôn là idx + 1)
			}
		}

		nodeCount = flat_nodes.size();
		cout << "DFS NodeCount: " << nodeCount << "\n";
	}

	void flat_bvh_to_bvh4(node* root) {
		if (!root) return;
		flat_nodes4.clear();

		// Hàng đợi lưu các nút BVH2 gốc sẽ trở thành các nút BVH4 lớn
		// và chỉ số vị trí tương ứng của chúng trong mảng phẳng `flat_nodes4` để cập nhật index sau
		struct QueueEntry {
			node* bvh2_node;
			int parent_flat4_idx;
			int child_slot; // 0, 1, 2, hoặc 3
		};

		queue<QueueEntry> q;
		q.push({ root, -1, 0 });

		while (!q.empty()) {
			auto current_entry = q.front();
			q.pop();

			node* src_node = current_entry.bvh2_node;

			// Nếu nút gốc BVH2 vốn đã là lá, ta tạo 1 nút BVH4 có 1 con là lá, 3 con còn lại trống
			// Tuy nhiên trường hợp phổ biến là nút TRONG. Ta cần nhặt tối đa 4 nút con từ cây BVH2:
			node* children2[4] = { nullptr, nullptr, nullptr, nullptr };
			int child_count = 0;

			if (src_node->leaf) {
				children2[0] = src_node;
				child_count = 1;
			}
			else {
				// Thuật toán thu thập 4 nút con hợp lệ từ 2 tầng của cây BVH2 hiện tại
				// Tầng 1:
				node* c1 = src_node->nodes[0];
				node* c2 = src_node->nodes[1];

				// Đánh giá mở rộng sang Tầng 2 để lấy đủ 4 con
				if (c1 && !c1->leaf && c2 && !c2->leaf) {
					children2[0] = c1->nodes[0];
					children2[1] = c1->nodes[1];
					children2[2] = c2->nodes[0];
					children2[3] = c2->nodes[1];
					child_count = 4;
				}
				else if (c1 && !c1->leaf) {
					children2[0] = c1->nodes[0];
					children2[1] = c1->nodes[1];
					children2[2] = c2;
					child_count = 3;
				}
				else if (c2 && !c2->leaf) {
					children2[0] = c1;
					children2[1] = c2->nodes[0];
					children2[2] = c2->nodes[1];
					child_count = 3;
				}
				else {
					children2[0] = c1;
					children2[1] = c2;
					child_count = 2;
				}
			}

			// Tạo nút BVH4 mới
			FlatNode4 n4;
			// Khởi tạo mặc định cho các kênh trống
			n4.bmin_x = n4.bmin_y = n4.bmin_z = vec4(0.0f);
			n4.bmax_x = n4.bmax_y = n4.bmax_z = vec4(0.0f);
			n4.child_or_start = ivec4(-1);
			n4.triangle_count = ivec4(-1); // -1 nghĩa là slot trống hoàn toàn

			int my_flat4_idx = flat_nodes4.size();
			flat_nodes4.push_back(n4);

			// Cập nhật con trỏ liên kết từ nút cha BVH4 trước đó tới nút hiện tại
			if (current_entry.parent_flat4_idx != -1) {
				flat_nodes4[current_entry.parent_flat4_idx].child_or_start[current_entry.child_slot] = my_flat4_idx;
			}

			// Đổ dữ liệu từ các nút BVH2 thu thập được vào 4 slot của nút BVH4
			for (int i = 0; i < child_count; ++i) {
				node* c = children2[i];
				if (!c) continue;

				// Nạp tọa độ Bound của từng hộp con vào các kênh tương ứng (x->[0], y->[1], z->[2], w->[3])
				flat_nodes4[my_flat4_idx].bmin_x[i] = c->box.bbox[0].x;
				flat_nodes4[my_flat4_idx].bmin_y[i] = c->box.bbox[0].y;
				flat_nodes4[my_flat4_idx].bmin_z[i] = c->box.bbox[0].z;

				flat_nodes4[my_flat4_idx].bmax_x[i] = c->box.bbox[1].x;
				flat_nodes4[my_flat4_idx].bmax_y[i] = c->box.bbox[1].y;
				flat_nodes4[my_flat4_idx].bmax_z[i] = c->box.bbox[1].z;

				if (c->leaf) {
					flat_nodes4[my_flat4_idx].child_or_start[i] = c->start;
					flat_nodes4[my_flat4_idx].triangle_count[i] = c->range; // > 0 là nút lá
				}
				else {
					flat_nodes4[my_flat4_idx].triangle_count[i] = 0; // 0 là nút trong, sẽ cập nhật child_index sau khi duyệt tới nó
					// Đẩy nút con này vào hàng đợi để tiếp tục bẻ thành nút BVH4 tiếp theo
					q.push({ c, my_flat4_idx, i });
				}
			}
		}
		nodeCount = flat_nodes4.size();
		cout << "Gộp BVH4 hoàn tất. Tổng số nút BVH4: " << flat_nodes4.size() << "\n";
	}

	void flat_bvh_to_bvh4_dfs(node* root) {
		if (!root) return;
		flat_nodes4.clear();

		// Cấu trúc Stack lưu thông tin duyệt DFS
		struct StackEntry {
			node* bvh2_node;
			int parent_flat4_idx; // Vị trí nút cha trong flat_nodes4 để quay lại sửa index
			int child_slot;        // Slot con (0, 1, 2, hoặc 3) của nút cha
		};

		vector<StackEntry> stk;
		stk.push_back({ root, -1, 0 });

		while (!stk.empty()) {
			auto current_entry = stk.back();
			stk.pop_back();

			node* src_node = current_entry.bvh2_node;

			// NẾU LÀ NÚT LÁ: Không cần tạo nút BVH4 bọc ngoài nữa, 
			// vì cấu trúc BVH4 quản lý lá trực tiếp từ slot của nút cha.
			if (src_node->leaf) {
				if (current_entry.parent_flat4_idx != -1) {
					int p_idx = current_entry.parent_flat4_idx;
					int slot = current_entry.child_slot;
					flat_nodes4[p_idx].child_or_start[slot] = src_node->start;
					flat_nodes4[p_idx].triangle_count[slot] = src_node->range; // > 0 là nút lá
				}
				continue;
			}

			// NÚT TRONG: Thu thập tối đa 4 nút con từ 2 tầng của cây BVH2
			node* children2[4] = { nullptr, nullptr, nullptr, nullptr };
			int child_count = 0;

			node* c1 = src_node->nodes[0];
			node* c2 = src_node->nodes[1];

			if (c1 && !c1->leaf && c2 && !c2->leaf) {
				children2[0] = c1->nodes[0];
				children2[1] = c1->nodes[1];
				children2[2] = c2->nodes[0];
				children2[3] = c2->nodes[1];
				child_count = 4;
			}
			else if (c1 && !c1->leaf) {
				children2[0] = c1->nodes[0];
				children2[1] = c1->nodes[1];
				children2[2] = c2;
				child_count = 3;
			}
			else if (c2 && !c2->leaf) {
				children2[0] = c1;
				children2[1] = c2->nodes[0];
				children2[2] = c2->nodes[1];
				child_count = 3;
			}
			else {
				children2[0] = c1;
				children2[1] = c2;
				child_count = 2;
			}

			// Tạo nút BVH4 mới
			FlatNode4 n4;
			n4.bmin_x = n4.bmin_y = n4.bmin_z = vec4(0.0f);
			n4.bmax_x = n4.bmax_y = n4.bmax_z = vec4(0.0f);
			n4.child_or_start = ivec4(-1);
			n4.triangle_count = ivec4(-1); // -1 nghĩa là slot trống

			int my_flat4_idx = flat_nodes4.size();
			flat_nodes4.push_back(n4);

			// CẬP NHẬT LIÊN KẾT TỪ NÚT CHA (BACKPATCHING CHÍNH XÁC)
			if (current_entry.parent_flat4_idx != -1) {
				int p_idx = current_entry.parent_flat4_idx;
				int slot = current_entry.child_slot;
				// Gán child_or_start trỏ tới vị trí của nút BVH4 hiện tại trong mảng phẳng
				flat_nodes4[p_idx].child_or_start[slot] = my_flat4_idx;
				flat_nodes4[p_idx].triangle_count[slot] = 0; // 0 nghĩa là nút trong
			}

			// Đổ dữ liệu từ BVH2 vào nút BVH4 mới và đẩy tiếp vào Stack theo thứ tự ngược (để pop ra đúng DFS)
			for (int i = child_count - 1; i >= 0; --i) {
				node* c = children2[i];
				if (!c) continue;

				flat_nodes4[my_flat4_idx].bmin_x[i] = c->box.bbox[0].x;
				flat_nodes4[my_flat4_idx].bmin_y[i] = c->box.bbox[0].y;
				flat_nodes4[my_flat4_idx].bmin_z[i] = c->box.bbox[0].z;

				flat_nodes4[my_flat4_idx].bmax_x[i] = c->box.bbox[1].x;
				flat_nodes4[my_flat4_idx].bmax_y[i] = c->box.bbox[1].y;
				flat_nodes4[my_flat4_idx].bmax_z[i] = c->box.bbox[1].z;

				if (c->leaf) {
					// Nếu con lượm được là nút lá, điền trực tiếp dữ liệu tam giác vào slot luôn
					flat_nodes4[my_flat4_idx].child_or_start[i] = c->start;
					flat_nodes4[my_flat4_idx].triangle_count[i] = c->range;
				}
				else {
					// Nếu là nút trong, đẩy vào stack để xử lý ở các bước tiếp theo
					stk.push_back({ c, my_flat4_idx, i });
				}
			}
		}

		nodeCount = flat_nodes4.size();
		cout << "Gop BVH4 (DFS) hoan tat. Tong so nut: " << nodeCount << "\n";
	}

	int count_nodes(node* n) {
		if (!n) return 0;
		if (n->leaf) return 1;
		return 1 + count_nodes(n->nodes[0]) + count_nodes(n->nodes[1]);
	}

	void flat_bvh_dfs_stackless2(node* n)
	{
		if (!n) return;

		flat_nodes.clear();
		flat_nodes.resize(count_nodes(n)); // C?p phát c?ng tr??c dung l??ng ?? tránh reallocate vùng nh? c?a vector

		// S? d?ng c? ch? DFS Backpatching chính xác tuy?t ??i
		struct StackEntry {
			node* ptr;
			int skip_idx;
		};

		std::vector<StackEntry> stk;
		stk.push_back({ n, -1 });

		int flat_idx_counter = 0;

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int my_idx = flat_idx_counter++;

			FlatNode flatnode;
			if (curr->leaf)
			{
				// Nút lá: box_min.w l?u v? trí start c?a tam giác, box_max.w l?u range (>0)
				flatnode.box_min = vec4(curr->box.bbox[0], (float)curr->start);
				flatnode.box_max = vec4(curr->box.bbox[1], (float)curr->range);
			}
			else
			{
				// Nút trong: box_min.w l?u skip_index cho GPU nh?y lúc tr??t h?p, box_max.w = 0
				flatnode.box_min = vec4(curr->box.bbox[0], (float)top.skip_idx);
				flatnode.box_max = vec4(curr->box.bbox[1], 0.0f);
			}

			flat_nodes[my_idx] = flatnode;

			if (!curr->leaf)
			{
				// L?y kích th??c cây con bên trái ?? bi?t chính xác v? trí con ph?i n?m ? ?âu trong m?ng ph?ng
				int left_subtree_size = count_nodes(curr->nodes[0]);
				int right_child_flat_idx = my_idx + 1 + left_subtree_size;

				// ??y con ph?i vào tr??c (N?u con ph?i tr??t, m?c tiêu c?a nó là skip_idx c?a cha)
				if (curr->nodes[1]) stk.push_back({ curr->nodes[1], top.skip_idx });

				// ??y con trái vào sau (N?u con trái tr??t, m?c tiêu c?a nó là ??u c?a nhánh con ph?i!)
				if (curr->nodes[0]) stk.push_back({ curr->nodes[0], right_child_flat_idx });
			}
		}

		nodeCount = flat_nodes.size();
		std::cout << "Stackless SBVH NodeCount: " << nodeCount << "\n";
	}

	void flat_bvh_dfs_stackless(node* n)
	{
		if (!n) return;

		flat_nodes.clear();
		//flat_nodes.resize(count_nodes(n));

		struct StackEntry {
			node* ptr;
			int skip_idx; // Ch? s? s? nh?y t?i n?u nút này ho?c c?m con c?a nó b? tr??t
		};

		vector<StackEntry> stk;
		// Nút g?c (Root) n?u tr??t thì thoát vòng l?p, nên skip_idx = -1
		stk.push_back({ n, -1 });

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			int current_idx = flat_nodes.size();

			FlatNode flatnode;
			if (curr->leaf)
			{
				// Nút lá: box_min.w l?u skip_index, box_max.w l?u range s? l??ng tam giác
				flatnode.box_min = vec4(curr->box.bbox[0], (float)top.skip_idx);
				flatnode.box_max = vec4(curr->box.bbox[1], (float)curr->range);
			}
			else
			{
				// Nút trong: box_min.w t?m th?i ?? 0, box_max.w = 0 (?ánh d?u nút trong)
				flatnode.box_min = vec4(curr->box.bbox[0], (float)top.skip_idx);
				flatnode.box_max = vec4(curr->box.bbox[1], 0.0f);
			}

			flat_nodes.emplace_back(flatnode);

			if (!curr->leaf)
			{
				// V?i c?u trúc Stackless DFS:
				// - N?u CON PH?I tr??t, nó s? nh?y t?i m?c tiêu 'skip_idx' c?a nút CHA.
				// - N?u CON TRÁI tr??t, nó s? nh?y th?ng sang ??u c?a nhánh CON PH?I.
				// Vì ta ch?a bi?t index chính xác c?a con ph?i (do ch?a push vào m?ng ph?ng),
				// ta dùng m?t m?o tính toán d?a trên kích th??c ??c l??ng ho?c qu?n lý liên k?t ng??c:

				// ?? ??n gi?n và chính xác tuy?t ??i mà không c?n ?? quy quay l?i s?a, 
				// Jacco Bikker g?i ý duy?t lu?ng stack l?u v? trí qu?n lý nh? sau:

				int right_child_skip_target = top.skip_idx;

				// ??y con ph?i vào tr??c
				stk.push_back({ curr->nodes[1], right_child_skip_target });

				// Con trái ???c ??y vào sau (x? lý tr??c). N?u con trái tr??t, 
				// m?c tiêu skip c?a nó chính là v? trí ID c?a con ph?i k? ti?p!
				// Trong m?ng ph?ng DFS, con ph?i luôn b?t ??u ngay sau khi toàn b? nhánh con trái k?t thúc.
				// Do ?ó ta c?n m?t ki?n trúc ?ánh d?u ch? s? ?? ?i?n ng??c (Backpatching).
			}
		}

		// --- TH?C HI?N BACKPATCHING TOÀN DI?N CHO SKIP POINTERS ---
		// Duy?t l?i c?u trúc cây ?? ?i?n chính xác ch? s? c?a con ph?i vào skip_index c?a con trái
		std::stack<std::pair<node*, int>> build_stk;
		build_stk.push({ n, -1 });
		int flat_idx_counter = 0;

		while (!build_stk.empty()) {
			auto top = build_stk.top();
			build_stk.pop();

			node* curr = top.first;
			int skip_target = top.second;

			flat_nodes[flat_idx_counter].box_min.w = (float)skip_target;
			int my_idx = flat_idx_counter;
			flat_idx_counter++;

			if (!curr->leaf) {
				// Tính toán v? trí ph?ng c?a con ph?i d?a trên kích th??c cây con bên trái
				int left_subtree_size = count_nodes(curr->nodes[0]);
				int right_child_flat_idx = my_idx + 1 + left_subtree_size;

				// Con ph?i tr??t -> nh?y t?i m?c tiêu cha
				build_stk.push({ curr->nodes[1], skip_target });
				// Con trái tr??t -> nh?y sang ??u con ph?i
				build_stk.push({ curr->nodes[0], right_child_flat_idx });
			}
		}

		nodeCount = flat_nodes.size();
		cout << "Stackless SBVH NodeCount: " << nodeCount << "\n";
	}

	void flat_bvh_dfs_parent_pointers(node* n)
	{
		if (!n) return;

		flat_nodes.clear();
		// ??m và c?p phát ?? b? nh? tr??c
		int total_nodes = count_nodes(n);
		flat_nodes.resize(total_nodes);

		struct StackEntry {
			node* ptr;
			int parent_idx;
			int is_right_child;
		};

		std::vector<StackEntry> stk;
		// Nút g?c không có cha, gán parent = -1
		stk.push_back({ n, -1, 0 });

		int flat_idx_counter = 0;

		// B? nh? ??m t?m th?i ?? d?ng cây
		struct TempNode {
			node* ptr;
			int my_flat_idx;
		};
		std::vector<TempNode> process_order;

		// Duy?t l??t 1: ??y các nút vào m?ng ph?ng theo chu?n c?p con n?m c?nh nhau
		// ?? ??t hi?u n?ng cao nh?t, ta duy?t theo hàng ??i ho?c dùng c?u trúc s?p x?p c?a Jacco Bikker
		// D??i ?ây là cách gán ch? s? Parent chu?n xác tuy?t ??i:

		// T?o c?u trúc cây DFS ph?ng chu?n hóa c?a b?n
		flat_bvh_dfs4(n); // Ch?y hàm ph?ng c? c?a b?n ?? l?y m?ng flat_nodes thô

		// Duy?t l??t 2: ?i?n con tr? cha (Parent Pointers) d?a trên liên k?t index
		std::stack<std::pair<node*, int>> parent_stk;
		parent_stk.push({ n, -1 }); // {node, parent_flat_index}
		int current_idx = 0;

		while (!parent_stk.empty()) {
			auto top = parent_stk.top();
			parent_stk.pop();

			node* curr = top.first;
			int p_idx = top.second;

			// Ghi ch? s? cha vào kênh .w c?a box_min
			flat_nodes[current_idx].box_min.w = (float)p_idx;
			int my_idx = current_idx;
			current_idx++;

			if (!curr->leaf) {
				// Cú pháp DFS c?a b?n: Con trái luôn là my_idx + 1
				// Tính toán v? trí con ph?i d?a trên kích th??c cây con bên trái
				int left_size = count_nodes(curr->nodes[0]);
				int right_flat_idx = my_idx + 1 + left_size;

				// ??y con ph?i vào tr??c, con trái vào sau ?? con trái ???c x? lý tr??c
				parent_stk.push({ curr->nodes[1], my_idx });
				parent_stk.push({ curr->nodes[0], my_idx });
			}
		}

		nodeCount = flat_nodes.size();
		std::cout << "Parent-Pointer SBVH NodeCount: " << nodeCount << "\n";
	}

	void flat_bvh_adjacent_children(node* n)
	{
		if (!n) return;

		flat_nodes.clear();

		// Kh?i t?o nút g?c t?i v? trí 0
		FlatNode root_flat;
		if (n->leaf) {
			root_flat.box_min = vec4(n->box.bbox[0], n->start);
			root_flat.box_max = vec4(n->box.bbox[1], n->range);
		}
		else {
			root_flat.box_min = vec4(n->box.bbox[0], 0); // S? c?p nh?t left_child_idx vào .w sau
			root_flat.box_max = vec4(n->box.bbox[1], 0);
		}
		flat_nodes.push_back(root_flat);

		// Stack l?u các nút trong (internal node) c?n x? lý các con c?a nó
		// parent_flat_idx: v? trí c?a nút cha trong m?ng flat_nodes
		struct StackEntry {
			node* ptr;
			int parent_flat_idx;
		};

		vector<StackEntry> stk;
		if (!n->leaf) {
			stk.push_back({ n, 0 });
		}

		while (!stk.empty())
		{
			auto top = stk.back();
			stk.pop_back();

			node* curr = top.ptr;
			node* left_child = curr->nodes[0];
			node* right_child = curr->nodes[1];

			// V? trí b?t ??u c?a 2 nút con trong m?ng flat_nodes (luôn k? nhau)
			int left_child_idx = flat_nodes.size();

			// 1. C?p nh?t ch? s? con trái vào nút cha (Con ph?i t? hi?u là left_child_idx + 1)
			flat_nodes[top.parent_flat_idx].box_min.w = (float)left_child_idx;

			// 2. T?o và ??y ??ng th?i c? 2 nút con vào m?ng flat_nodes
			FlatNode flat_left, flat_right;

			// C?u hình cho con TRÁI
			if (left_child->leaf) {
				flat_left.box_min = vec4(left_child->box.bbox[0], left_child->start);
				flat_left.box_max = vec4(left_child->box.bbox[1], left_child->range);
			}
			else {
				flat_left.box_min = vec4(left_child->box.bbox[0], 0);
				flat_left.box_max = vec4(left_child->box.bbox[1], 0);
			}

			// C?u hình cho con PH?I
			if (right_child->leaf) {
				flat_right.box_min = vec4(right_child->box.bbox[0], right_child->start);
				flat_right.box_max = vec4(right_child->box.bbox[1], right_child->range);
			}
			else {
				flat_right.box_min = vec4(right_child->box.bbox[0], 0);
				flat_right.box_max = vec4(right_child->box.bbox[1], 0);
			}

			flat_nodes.push_back(flat_left);
			flat_nodes.push_back(flat_right);

			// 3. ??y các nút con vào stack ?? ti?p t?c duy?t sâu xu?ng d??i
			// ??y con PH?I vào tr??c, con TRÁI vào sau ?? khi pop thì x? lý bên TRÁI tr??c (gi? tính ch?t DFS)
			if (!right_child->leaf) {
				stk.push_back({ right_child, left_child_idx + 1 });
			}
			if (!left_child->leaf) {
				stk.push_back({ left_child, left_child_idx });
			}
		}

		nodeCount = flat_nodes.size();
		cout << "Adjacent BVH NodeCount: " << nodeCount << "\n";
	}

	void convert_to_soa()
	{
		for (auto& v : flat_nodes)
		{
			vec4 bmin = v.box_min;
			vec4 bmax = v.box_max;

			bvh_soa.min_x.emplace_back(bmin.x);
			bvh_soa.min_y.emplace_back(bmin.y);
			bvh_soa.min_z.emplace_back(bmin.z);

			bvh_soa.start.emplace_back(bmin.w);

			bvh_soa.max_x.emplace_back(bmax.x);
			bvh_soa.max_y.emplace_back(bmax.y);
			bvh_soa.max_z.emplace_back(bmax.z);

			bvh_soa.range.emplace_back(bmax.w);
		}

		//cout << "min_x " << bvh_soa.min_x.size() << "\n";

		//vector<FlatNode>().swap(flat_nodes);
	}

	void convert_to_center_form()
	{
		for (auto& v : flat_nodes)
		{
			vec3 box_min_3 = vec3(v.box_min.x, v.box_min.y, v.box_min.z);
			vec3 box_max_3 = vec3(v.box_max.x, v.box_max.y, v.box_max.z);

			vec3 center = (box_min_3 + box_max_3) * 0.5f;
			vec3 extend = (box_max_3 - box_min_3) * 0.5f;

			v.box_min.x = center.x;
			v.box_min.y = center.y;
			v.box_min.z = center.z;

			v.box_max.x = extend.x;
			v.box_max.y = extend.y;
			v.box_max.z = extend.z;
		}
	}

	void clean_up()
	{
		vector<Reference>().swap(references);

		vector<int>().swap(triangle_indices);

		vector<BBox>().swap(m_rightBound);

		vector<FlatNode>().swap(flat_nodes);

		vector<FlatNode4>().swap(flat_nodes4);
	}
	void delete_bvh(node*& n)
	{
		if (n == NULL)
			return;
		delete_bvh(n->nodes[0]);
		delete_bvh(n->nodes[1]);

		delete(n);
	}

	int height(node*& n)
	{
		if (n == NULL)
			return -1;

		int left_height = height(n->nodes[0]);
		int right_height = height(n->nodes[1]);

		return maxf(left_height, right_height) + 1;
	}

};

#endif // !_SBVH2_H_

