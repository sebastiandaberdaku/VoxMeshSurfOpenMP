#ifndef MESH_H
#define MESH_H

#include "primitives.h"

#include <iostream>
using std::cout;
using std::endl;

#include <fstream>
using std::ifstream;
using std::ofstream;

#include <iomanip>
using std::setiosflags;

#include <ios>
using std::ios_base;
using std::ios;

#include <set>
using std::set;

#include <vector>
using std::vector;

#include <limits>
using std::numeric_limits;

#include <cstring> // for memcpy()
#include <cctype>


class mesh
{
public:
	vector<vertex_3> vertices;
	vector<indexed_triangle> triangles;
	vector< vector<size_t> > vertex_to_vertex_indices;
	vector< vector<size_t> > vertex_to_triangle_indices;
	vector<vertex_3> vertex_normals;
	vector<vertex_3> triangle_normals;

	void clear(void);
	bool operator==(const mesh &right);
	bool operator!=(const mesh &right);

	void fix_cracks(void);

	bool load_from_binary_stereo_lithography_file(const char *const file_name, const bool generate_normals = true, const size_t buffer_width = 65536);
	bool save_to_binary_stereo_lithography_file(const char *const file_name, const size_t buffer_width = 65536);
	bool save_to_povray_mesh2_file(const char *const file_name, const bool write_vertex_normals = false);

	float get_max_extent(void);
	float set_max_extent(const float target_max_extent);
	float get_triangle_area(const size_t tri_index);
	float get_vertex_neighbourhood_area(const size_t vertex_index);
	float get_area(void);
	float set_area(const float target_area);
	float get_triangle_volume(const size_t tri_index);
	float get_volume(void);
	float set_volume(const float target_volume);

	// See: Geometric Signal Processing on Polygonal Meshes by G. Taubin
	void laplace_smooth(const float scale);
	void laplace_smooth_fujiwara(const float scale);
	void laplace_smooth_cn(const float scale);
	void taubin_smooth(const float lambda, const float mu, const size_t steps);
	void taubin_smooth_fujiwara(const float lambda, const float mu, const size_t steps);
	void taubin_smooth_cn(const float lambda, const float mu, const size_t steps);

protected:
	void generate_vertex_normals(void);
	void generate_triangle_normals(void);
	void generate_vertex_and_triangle_normals(void);
	void regenerate_vertex_and_triangle_normals_if_exists(void);
	template<typename T> void eliminate_vector_duplicates(vector<T> &v);
	bool merge_vertex_pair(const size_t keeper, const size_t goner);
};


#endif
