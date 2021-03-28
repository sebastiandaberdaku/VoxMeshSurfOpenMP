#include "mesh.h"

void mesh::clear(void)
{
	triangles.clear();
	vertices.clear();
	vertex_to_triangle_indices.clear();
	vertex_to_vertex_indices.clear();
	vertex_normals.clear();
	triangle_normals.clear();
}

bool mesh::operator==(const mesh &right)
{
	if(triangles == right.triangles &&	vertices == right.vertices && vertex_to_triangle_indices == right.vertex_to_triangle_indices &&	vertex_to_vertex_indices == right.vertex_to_vertex_indices)
		return true;

	return false;
}

bool mesh::operator!=(const mesh &right)
{
	return !(*this == right);
}

void mesh::fix_cracks(void)
{
	cout << "Finding cracks" << endl;

	// Find edges that belong to only one triangle.
	set<ordered_indexed_edge> problem_edges_set;
	size_t problem_edge_id = 0;

	// For each vertex.
	for(size_t i = 0; i < vertices.size(); i++)
	{
		// For each edge.
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t triangle_count = 0;
			size_t neighbour_j = vertex_to_vertex_indices[i][j];

			// Find out which triangles are shared by this edge.
			for(size_t k = 0; k < vertex_to_triangle_indices[i].size(); k++)
			{
				for(size_t l = 0; l < vertex_to_triangle_indices[neighbour_j].size(); l++)
				{
					size_t tri0_index = vertex_to_triangle_indices[i][k];
					size_t tri1_index = vertex_to_triangle_indices[neighbour_j][l];

					if(tri0_index == tri1_index)
					{
						triangle_count++;
						break;
					}
				}
			}

			// Found a problem edge.
			if(2 != triangle_count)
			{
				indexed_vertex_3 v0;
				v0.index = i;
				v0.x = vertices[i].x;
				v0.y = vertices[i].y;
				v0.z = vertices[i].z;

				indexed_vertex_3 v1;
				v1.index = neighbour_j;
				v1.x = vertices[neighbour_j].x;
				v1.y = vertices[neighbour_j].y;
				v1.z = vertices[neighbour_j].z;

				ordered_indexed_edge problem_edge(v0, v1);

				if(problem_edges_set.end() == problem_edges_set.find(problem_edge))
				{
					problem_edge.id = problem_edge_id++;
					problem_edges_set.insert(problem_edge);
				}
			} // End of: Found a problem edge.
		} // End of: For each edge.
	} // End of: For each vertex.

	if(0 == problem_edges_set.size())
	{
		cout << "No cracks found -- the mesh seems to be in good condition" << endl;
		return;
	}

	cout << "Found " << problem_edges_set.size() << " problem edges" << endl;

	if(0 != problem_edges_set.size() % 2)
	{
		cout << "Error -- the number of problem edges must be an even number (perhaps the mesh has holes?). Aborting." << endl;
		return;
	}

	// Make a copy of the set into a vector because the edge matching will
	// run a bit faster while looping through a vector by index vs looping through
	// a set by iterator.
	vector<ordered_indexed_edge> problem_edges_vec(problem_edges_set.begin(), problem_edges_set.end());
	vector<bool> processed_problem_edges(problem_edges_set.size(), false);
	problem_edges_set.clear();

	set<ordered_size_t_pair> merge_vertices;

	cout << "Pairing problem edges" << endl;

	// Each problem edge is practically a duplicate of some other, but not quite exactly.
	// So, find the closest match for each problem edge.
	for(size_t i = 0; i < problem_edges_vec.size(); i++)
	{
		// This edge has already been matched up, so skip it.
		if(true == processed_problem_edges[problem_edges_vec[i].id])
			continue;

		float closest_dist_sq = numeric_limits<float>::max();
		size_t closest_problem_edges_vec_index = 0;

		for(size_t j = i + 1; j < problem_edges_vec.size(); j++)
		{
			// Note: Don't check to see if this edge has already been matched up.
			// Doing so will actually only slow this down further. Perhaps vector<bool> is a bit of a sloth?

			float dist_sq = problem_edges_vec[i].centre_point.distance_sq(problem_edges_vec[j].centre_point);

			if(dist_sq < closest_dist_sq)
			{
				closest_dist_sq = dist_sq;
				closest_problem_edges_vec_index = j;
			}
		}

		processed_problem_edges[problem_edges_vec[i].id] = true;
		processed_problem_edges[problem_edges_vec[closest_problem_edges_vec_index].id] = true;

		// If edge 0 vertex 0 is further in space from edge 1 vertex 0 than from edge 1 vertex 1,
		// then swap the indices on the edge 1 -- this makes sure that the edges are not pointing
		// in opposing directions.
		if(vertices[problem_edges_vec[i].indices[0]].distance_sq(vertices[problem_edges_vec[closest_problem_edges_vec_index].indices[0]]) > vertices[problem_edges_vec[i].indices[0]].distance_sq(vertices[problem_edges_vec[closest_problem_edges_vec_index].indices[1]]))
		{
			size_t temp = problem_edges_vec[closest_problem_edges_vec_index].indices[0];
			problem_edges_vec[closest_problem_edges_vec_index].indices[0] = problem_edges_vec[closest_problem_edges_vec_index].indices[1];
			problem_edges_vec[closest_problem_edges_vec_index].indices[1] = temp;
		}

		// If the first indices aren't already the same, then merge them.
		if(problem_edges_vec[i].indices[0] != problem_edges_vec[closest_problem_edges_vec_index].indices[0])
			merge_vertices.insert(ordered_size_t_pair(problem_edges_vec[i].indices[0], problem_edges_vec[closest_problem_edges_vec_index].indices[0]));

		// If the second indices aren't already the same, then merge them.
		if(problem_edges_vec[i].indices[1] != problem_edges_vec[closest_problem_edges_vec_index].indices[1])
			merge_vertices.insert(ordered_size_t_pair(problem_edges_vec[i].indices[1], problem_edges_vec[closest_problem_edges_vec_index].indices[1]));
	}

	cout << "Merging " << merge_vertices.size() << " vertex pairs" << endl;

	for(set<ordered_size_t_pair>::const_iterator ci = merge_vertices.begin(); ci != merge_vertices.end(); ci++)
		merge_vertex_pair(ci->indices[0], ci->indices[1]);

	// Recalculate normals, if necessary.
	regenerate_vertex_and_triangle_normals_if_exists();
}

bool mesh::load_from_binary_stereo_lithography_file(const char *const file_name, const bool generate_normals, const size_t buffer_width)
{
	clear();

	cout << "Reading file: " << file_name << endl;

	ifstream in(file_name, ios_base::binary);

	if(in.fail())
		return false;

	const size_t header_size = 80;
	vector<char> buffer(header_size, 0);
	unsigned int num_triangles = 0; // Must be 4-byte unsigned int.

	// Read header.
	in.read(reinterpret_cast<char *>(&(buffer[0])), header_size);

	if(header_size != in.gcount())
		return false;

	if( 's' == tolower(buffer[0]) &&
		'o' == tolower(buffer[1]) && 
		'l' == tolower(buffer[2]) && 
		'i' == tolower(buffer[3]) && 
		'd' == tolower(buffer[4]) )
	{
		cout << "Encountered ASCII STL file header -- aborting." << endl;
		return false;
	}

	// Read number of triangles.
	in.read(reinterpret_cast<char *>(&num_triangles), sizeof(unsigned int));

	if(sizeof(unsigned int) != in.gcount())
		return false;

	triangles.resize(num_triangles);

	cout << "Triangles:    " << triangles.size() << endl;

	// Enough bytes for twelve 4-byte floats plus one 2-byte integer, per triangle.
	const size_t per_triangle_data_size = (12*sizeof(float) + sizeof(short unsigned int));
	const size_t buffer_size = per_triangle_data_size * buffer_width;
	buffer.resize(buffer_size, 0);

	size_t num_triangles_remaining = triangles.size();
	size_t tri_index = 0;
	set<indexed_vertex_3> vertex_set;

	while(num_triangles_remaining > 0)
	{
		size_t num_triangles_to_read = buffer_width;

		if(num_triangles_remaining < buffer_width)
			num_triangles_to_read = num_triangles_remaining;

		size_t data_size = per_triangle_data_size*num_triangles_to_read;

		in.read(reinterpret_cast<char *>(&buffer[0]), data_size);

		if(data_size != in.gcount() || in.fail())
			return false;

		num_triangles_remaining -= num_triangles_to_read;

		// Use a pointer to assist with the copying.
		// Should probably use std::copy() instead, but memcpy() does the trick, so whatever...
		char *cp = &buffer[0];

		for(size_t i = 0; i < num_triangles_to_read; i++)
		{
			// Skip face normal. We will calculate them manually later.
			cp += 3*sizeof(float);

			// For each of the three vertices in the triangle.
			for(short unsigned int j = 0; j < 3; j++)
			{
				indexed_vertex_3 v;

				// Get vertex components.
				memcpy(&v.x, cp, sizeof(float)); cp += sizeof(float);
				memcpy(&v.y, cp, sizeof(float)); cp += sizeof(float);
				memcpy(&v.z, cp, sizeof(float)); cp += sizeof(float);

				// Look for vertex in set.
				set<indexed_vertex_3>::const_iterator find_iter = vertex_set.find(v);

				// If vertex not found in set...
				if(vertex_set.end() == find_iter)
				{
					// Assign new vertices index
					v.index = vertices.size();

					// Add vertex to set
					vertex_set.insert(v);

					// Add vertex to vector
					vertex_3 indexless_vertex;
					indexless_vertex.x = v.x;
					indexless_vertex.y = v.y;
					indexless_vertex.z = v.z;
					vertices.push_back(indexless_vertex);

					// Assign vertex index to triangle
					triangles[tri_index].vertex_indices[j] = v.index;

					// Add triangle index to vertex
					vector<size_t> tri_indices;
					tri_indices.push_back(tri_index);
					vertex_to_triangle_indices.push_back(tri_indices);
				}
				else
				{
					// Assign existing vertex index to triangle
					triangles[tri_index].vertex_indices[j] = find_iter->index;

					// Add triangle index to vertex
					vertex_to_triangle_indices[find_iter->index].push_back(tri_index);
				}
			}

			// Skip attribute.
			cp += sizeof(short unsigned int);

			tri_index++;
		}
	}

	vertex_to_vertex_indices.resize(vertices.size());

	for(size_t i = 0; i < vertex_to_triangle_indices.size(); i++)
	{
		// Use a temporary set to avoid duplicates.
		set<size_t> vertex_to_vertex_indices_set;

		for(size_t j = 0; j < vertex_to_triangle_indices[i].size(); j++)
		{
			size_t tri_index = vertex_to_triangle_indices[i][j];

			for(size_t k = 0; k < 3; k++)
				if(i != triangles[tri_index].vertex_indices[k]) // Don't add current vertex index to its own adjacency list.
					vertex_to_vertex_indices_set.insert(triangles[tri_index].vertex_indices[k]);
		}

		// Copy to final vector.
		for(set<size_t>::const_iterator ci = vertex_to_vertex_indices_set.begin(); ci != vertex_to_vertex_indices_set.end(); ci++)
			vertex_to_vertex_indices[i].push_back(*ci);
	}

	cout << "Vertices:     " << triangles.size()*3 << " (of which " << vertices.size() << " are unique)" << endl;

	in.close();

	if(true == generate_normals)
	{
		cout << "Generating normals" << endl;
		generate_vertex_and_triangle_normals();
	}

	return true;
} 

bool mesh::save_to_binary_stereo_lithography_file(const char *const file_name, const size_t buffer_width)
{
	cout << "Writing file: " << file_name << endl;
	cout << "Triangles:    " << triangles.size() << endl;
	cout << "Vertices:     " << triangles.size()*3 << endl;

	if(0 == triangles.size())
		return false;

	// Write to file.
	ofstream out(file_name, ios_base::binary);

	if(out.fail())
		return false;

	const size_t header_size = 80;
	vector<char> buffer(header_size, 0);
	const unsigned int num_triangles = triangles.size(); // Must be 4-byte unsigned int.
	vertex_3 normal;

	// Write blank header.
	out.write(reinterpret_cast<const char *>(&(buffer[0])), header_size);

	// Write number of triangles.
	out.write(reinterpret_cast<const char *>(&num_triangles), sizeof(unsigned int));

	// Enough bytes for twelve 4-byte floats plus one 2-byte integer, per triangle.
	const size_t per_triangle_data_size = (12*sizeof(float) + sizeof(short unsigned int));
	const size_t buffer_size = per_triangle_data_size * buffer_width;
	buffer.resize(buffer_size, 0);

	// Use a pointer to assist with the copying.
	// Should probably use std::copy() instead, but memcpy() does the trick, so whatever...
	char *cp = &buffer[0];
	size_t buffer_count = 0;

	cout << "Writing " << per_triangle_data_size*triangles.size() / 1048576 << " MB of data to disk" << endl;

	for(size_t i = 0; i < triangles.size(); i++)
	{
		// Copy face normal if it's been calculated, otherwise manually calculate it.
		if(triangle_normals.size() == triangles.size())
		{
			normal = triangle_normals[i];
		}
		else
		{
			vertex_3 v0 = vertices[triangles[i].vertex_indices[1]] - vertices[triangles[i].vertex_indices[0]];
			vertex_3 v1 = vertices[triangles[i].vertex_indices[2]] - vertices[triangles[i].vertex_indices[0]];
			normal = v0.cross(v1);
			normal.normalize();
		}

		memcpy(cp, &normal.x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &normal.z, sizeof(float)); cp += sizeof(float);

		memcpy(cp, &vertices[triangles[i].vertex_indices[0]].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[0]].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[0]].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[1]].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[1]].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[1]].z, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[2]].x, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[2]].y, sizeof(float)); cp += sizeof(float);
		memcpy(cp, &vertices[triangles[i].vertex_indices[2]].z, sizeof(float)); cp += sizeof(float);

		cp += sizeof(short unsigned int);

		buffer_count++;

		// If buffer is full, write triangles in buffer to disk.
		if(buffer_count == buffer_width)
		{
			out.write(reinterpret_cast<const char *>(&buffer[0]), buffer_size);

			if(out.fail())
				return false;

			buffer_count = 0;
			cp = &buffer[0];
		}
	}

	// Write any remaining triangles in buffer to disk.
	// This will occur whenever triangles.size() % buffer_width != 0
	// (ie. when triangle count is not a multiple of buffer_width, which should happen almost all of the time).
	if(buffer_count > 0)
	{
		out.write(reinterpret_cast<const char *>(&buffer[0]), per_triangle_data_size*buffer_count);

		if(out.fail())
			return false;
	}

	out.close();

	return true;
} 

bool mesh::save_to_povray_mesh2_file(const char *const file_name, const bool write_vertex_normals)
{
	cout << "Triangle count: " << triangles.size() << endl;

	if(0 == triangles.size())
		return false;

	if(true == write_vertex_normals && vertex_normals.size() != vertices.size())
		generate_vertex_normals();

	// Write to file.
	ofstream out(file_name);

	if(out.fail())
		return false;

	out << setiosflags(ios_base::fixed);

	cout << "Writing data to " << file_name << endl;

	// Bump up output precision to help keep very small triangles from becoming degenerate.
	//out << setprecision(18);

	// Note: Some of these vertices may be rogue vertices that aren't referenced by triangles;
	// this occurs after cracks have been fixed. Whatever.
	out << " vertex_vectors" << endl;
	out << " {" << endl;
	out << "  " << vertices.size() << ',' << endl;

	for(size_t i = 0; i < vertices.size() - 1; i++)
		out << "  <" << vertices[i].x << ',' << vertices[i].y << ',' << vertices[i].z << ">," << endl;

	out << "  <" << vertices[vertices.size() - 1].x << ',' << vertices[vertices.size() - 1].y << ',' << vertices[vertices.size() - 1].z << '>' << endl;
	out << " }" << endl;

	if(true == write_vertex_normals)
	{
		out << " normal_vectors" << endl;
		out << " {" << endl;
		out << "  " << vertex_normals.size() << ',' << endl;

		for(size_t i = 0; i < vertex_normals.size() - 1; i++)
			out << "  <" << vertex_normals[i].x << ',' << vertex_normals[i].y << ',' << vertex_normals[i].z << ">," << endl;

		out << "  <" << vertex_normals[vertex_normals.size() - 1].x << ',' << vertex_normals[vertex_normals.size() - 1].y << ',' << vertex_normals[vertex_normals.size() - 1].z << '>' << endl;
		out << " }" << endl;
	}

	out << " face_indices" << endl;
	out << " {" << endl;
	out << "  " << triangles.size() << ',' << endl;

	for(size_t i = 0; i < triangles.size() - 1; i++)
		out << "  <" << triangles[i].vertex_indices[0] << ',' << triangles[i].vertex_indices[1] << ',' << triangles[i].vertex_indices[2] << ">," << endl;

	out << "  <" << triangles[triangles.size() - 1].vertex_indices[0] << ',' << triangles[triangles.size() - 1].vertex_indices[1] << ',' << triangles[triangles.size() - 1].vertex_indices[2]<< ">" << endl;
	out << " }" << endl;

	out.close();

	return true;
}

float mesh::get_max_extent(void)
{
	float current_x_min = numeric_limits<float>::max();
	float current_y_min = numeric_limits<float>::max();
	float current_z_min = numeric_limits<float>::max();
	float current_x_max = numeric_limits<float>::min();
	float current_y_max = numeric_limits<float>::min();
	float current_z_max = numeric_limits<float>::min();

	for(size_t i = 0; i < vertices.size(); i++)
	{
		if(vertices[i].x < current_x_min)
			current_x_min = vertices[i].x;

		if(vertices[i].x > current_x_max)
			current_x_max = vertices[i].x;

		if(vertices[i].y < current_y_min)
			current_y_min = vertices[i].y;

		if(vertices[i].y > current_y_max)
			current_y_max = vertices[i].y;

		if(vertices[i].z < current_z_min)
			current_z_min = vertices[i].z;

		if(vertices[i].z > current_z_max)
			current_z_max = vertices[i].z;
	}

	float current_x_extent = fabsf(current_x_min - current_x_max);
	float current_y_extent = fabsf(current_y_min - current_y_max);
	float current_z_extent = fabsf(current_z_min - current_z_max);

	float current_max_extent = current_x_extent;

	if(current_y_extent > current_max_extent)
		current_max_extent = current_y_extent;

	if(current_z_extent > current_max_extent)
		current_max_extent = current_z_extent;

	return current_max_extent;
}

float mesh::set_max_extent(const float target_max_extent)
{
	float current_max_extent = get_max_extent();

	float scale_value = target_max_extent / current_max_extent;

	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] *= scale_value;

	return get_max_extent();
}

float mesh::get_triangle_area(const size_t tri_index)
{
	if(tri_index >= triangles.size())
		return 0;

	vertex_3 a = vertices[triangles[tri_index].vertex_indices[1]] - vertices[triangles[tri_index].vertex_indices[0]];
	vertex_3 b = vertices[triangles[tri_index].vertex_indices[2]] - vertices[triangles[tri_index].vertex_indices[0]];
	vertex_3 cross = a.cross(b);

	return 0.5f*cross.length();
}

float mesh::get_vertex_neighbourhood_area(const size_t vertex_index)
{
	if(vertex_index >= vertices.size())
		return 0;

	float total_area = 0;

	for(size_t i = 0; i < vertex_to_triangle_indices[vertex_index].size(); i++)
		total_area += get_triangle_area(vertex_to_triangle_indices[vertex_index][i]);

	return total_area;
}

float mesh::get_area(void)
{
	float total_area = 0;

	for(size_t i = 0; i < triangles.size(); i++)
		total_area += get_triangle_area(i);

	return total_area;
}

float mesh::set_area(const float target_area)
{
	float current_area = get_area();

	if(0 == current_area)
		return 0;

	float scale_factor = pow(target_area / current_area, 0.5f);

	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] *= scale_factor;

	return get_area();
}

float mesh::get_triangle_volume(const size_t tri_index)
{
	if(tri_index >= triangles.size())
		return 0;

	vertex_3 a = vertices[triangles[tri_index].vertex_indices[0]];
	vertex_3 b = vertices[triangles[tri_index].vertex_indices[1]];
	vertex_3 c = vertices[triangles[tri_index].vertex_indices[2]];

	return a.dot(b.cross(c)) / 6.0f;
}

float mesh::get_volume(void)
{
	float total_volume = 0;

	for(size_t i = 0; i < triangles.size(); i++)
		total_volume += get_triangle_volume(i);

	return total_volume;
}

float mesh::set_volume(const float target_volume)
{
	float current_volume = get_volume();

	if(0 == current_volume)
		return 0;

	float scale_factor = pow(target_volume / current_volume, 1.0f/3.0f);

	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] *= scale_factor;

	return get_volume();
}

// This produces results that are practically identical to Meshlab
void mesh::laplace_smooth(const float scale)
{
	vector<vertex_3> displacements(vertices.size(), vertex_3(0, 0, 0));

	// Get per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
	{
		// Skip rogue vertices (which were probably made rogue during a previous
		// attempt to fix mesh cracks).
		if(0 == vertex_to_vertex_indices[i].size())
			continue;

		const float weight = 1.0f / static_cast<float>(vertex_to_vertex_indices[i].size());

		// Sum the displacements.
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t neighbour_j = vertex_to_vertex_indices[i][j];
			displacements[i] += (vertices[neighbour_j] - vertices[i])*weight;
		}
	}

	// Apply per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] += displacements[i]*scale;
}

void mesh::laplace_smooth_fujiwara(const float scale)
{
	vector<vertex_3> displacements(vertices.size(), vertex_3(0, 0, 0));

	// Get per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
	{
		// Skip rogue vertices (which were probably made rogue during a previous
		// attempt to fix mesh cracks).
		if(0 == vertex_to_vertex_indices[i].size())
			continue;

		vector<float> weights(vertex_to_vertex_indices[i].size(), 0.0f);

		// Calculate Fujiwara weights based on edge lengths.
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t neighbour_j = vertex_to_vertex_indices[i][j];

			float edge_length = vertices[i].distance(vertices[neighbour_j]);

			if(0 == edge_length)
				edge_length = numeric_limits<float>::epsilon();

			weights[j] = 1.0f / edge_length;
		}

		// Normalize the weights so that they sum up to 1.
		float s = 0;

		for(size_t j = 0; j < weights.size(); j++)
			s += weights[j];

		if(0 == s)
			s = numeric_limits<float>::epsilon();

		for(size_t j = 0; j < weights.size(); j++)
			weights[j] /= s;

		// Sum the displacements.
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t neighbour_j = vertex_to_vertex_indices[i][j];
			displacements[i] += (vertices[neighbour_j] - vertices[i])*weights[j];
		}
	}

	// Apply per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] += displacements[i]*scale;
}

void mesh::laplace_smooth_cn(const float scale)
{
	vector<vertex_3> displacements(vertices.size(), vertex_3(0, 0, 0));

	// Get per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
	{
		if(0 == vertex_to_vertex_indices[i].size())
			continue;

		vector<float> weights(vertex_to_vertex_indices[i].size(), 0.0f);

		size_t angle_error = 0;

		// For each vertex pair (ie. each edge),
		// calculate weight based on the two opposing angles (ie. curvature normal scheme).
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t angle_count = 0;

			size_t neighbour_j = vertex_to_vertex_indices[i][j];

			// Find out which two triangles are shared by the edge.
			for(size_t k = 0; k < vertex_to_triangle_indices[i].size(); k++)
			{
				for(size_t l = 0; l < vertex_to_triangle_indices[neighbour_j].size(); l++)
				{
					size_t tri0_index = vertex_to_triangle_indices[i][k];
					size_t tri1_index = vertex_to_triangle_indices[neighbour_j][l];

					// This will occur twice per edge.
					if(tri0_index == tri1_index)
					{
						// Find the third vertex in this triangle (the vertex that doesn't belong to the edge).
						for(size_t m = 0; m < 3; m++)
						{
							// This will occur once per triangle.
							if(triangles[tri0_index].vertex_indices[m] != i && triangles[tri0_index].vertex_indices[m] != neighbour_j)
							{
								size_t opp_vert_index = triangles[tri0_index].vertex_indices[m];

								// Get the angle opposite of the edge.
								vertex_3 a = vertices[i] - vertices[opp_vert_index];
								vertex_3 b = vertices[neighbour_j] - vertices[opp_vert_index];
								a.normalize();
								b.normalize();

								float dot = a.dot(b);

								if(-1 > dot)
									dot = -1;
								else if(1 < dot)
									dot = 1;

								float angle = acosf(dot);

								// Curvature normal weighting.
								float slope = tanf(angle);

								if(0 == slope)
									slope = numeric_limits<float>::epsilon();

								// Note: Some weights will be negative, due to obtuse triangles.
								// You may wish to do weights[j] += fabsf(1.0f / slope); here.
								weights[j] += 1.0f / slope;

								angle_count++;

								break;
							}
						}

						// Since we found a triangle match, we can skip to the first vertex's next triangle.
						break;
					}
				}
			} // End of: Find out which two triangles are shared by the vertex pair.

			if(angle_count != 2)
				angle_error++;

		} // End of: For each vertex pair (ie. each edge).

		if(angle_error != 0)
		{
			cout << "Warning: Vertex " << i << " belongs to " << angle_error << " edges that do not belong to two triangles (" << vertex_to_vertex_indices[i].size() - angle_error << " edges were OK)." << endl;
			cout << "Your mesh probably has cracks or holes in it." << endl;
		}

		// Normalize the weights so that they sum up to 1.
		float s = 0;

		// Note: Some weights will be negative, due to obtuse triangles.
		// You may wish to do s += fabsf(weights[j]); here.
		for(size_t j = 0; j < weights.size(); j++)
			s += weights[j];

		if(0 == s)
			s = numeric_limits<float>::epsilon();

		for(size_t j = 0; j < weights.size(); j++)
			weights[j] /= s;

		// Sum the displacements.
		for(size_t j = 0; j < vertex_to_vertex_indices[i].size(); j++)
		{
			size_t neighbour_j = vertex_to_vertex_indices[i][j];

			displacements[i] += (vertices[neighbour_j] - vertices[i])*weights[j];
		}
	}

	// Todo: Find out why there are cases where displacement is much, much, much larger than all edge lengths put together.

	// Apply per-vertex displacement.
	for(size_t i = 0; i < vertices.size(); i++)
		vertices[i] += displacements[i]*scale;
}



void mesh::taubin_smooth(const float lambda, const float mu, const size_t steps)
{
	cout << "Smoothing mesh using Taubin lambda|mu algorithm ";
	cout << "(inverse neighbour count weighting)" << endl;

	for(size_t s = 0; s < steps; s++)
	{
		cout << "Step " << s + 1 << " of " << steps << endl;

		laplace_smooth(lambda);
		laplace_smooth(mu);
	}

	// Recalculate normals, if necessary.
	regenerate_vertex_and_triangle_normals_if_exists();
}

void mesh::taubin_smooth_fujiwara(const float lambda, const float mu, const size_t steps)
{
	cout << "Smoothing mesh using Taubin lambda|mu algorithm ";
	cout << "(Fujiwara weighting)" << endl;

	for(size_t s = 0; s < steps; s++)
	{
		cout << "Step " << s + 1 << " of " << steps << endl;

		laplace_smooth_fujiwara(lambda);
		laplace_smooth_fujiwara(mu);
	}

	// Recalculate normals, if necessary.
	regenerate_vertex_and_triangle_normals_if_exists();
}

void mesh::taubin_smooth_cn(const float lambda, const float mu, const size_t steps)
{
	cout << "Smoothing mesh using Taubin lambda|mu algorithm ";
	cout << "(curvature normal weighting)" << endl;

	for(size_t s = 0; s < steps; s++)
	{
		cout << "Step " << s + 1 << " of " << steps << endl;

		laplace_smooth_cn(lambda);
		laplace_smooth_cn(mu);
	}

	// Recalculate normals, if necessary.
	regenerate_vertex_and_triangle_normals_if_exists();
}

void mesh::generate_vertex_normals(void)
{
	if(triangles.size() == 0 || vertices.size() == 0)
		return;

	vertex_normals.clear();
	vertex_normals.resize(vertices.size());

	for(size_t i = 0; i < triangles.size(); i++)
	{
		vertex_3 v0 = vertices[triangles[i].vertex_indices[1]] - vertices[triangles[i].vertex_indices[0]];
		vertex_3 v1 = vertices[triangles[i].vertex_indices[2]] - vertices[triangles[i].vertex_indices[0]];
		vertex_3 v2 = v0.cross(v1);

		vertex_normals[triangles[i].vertex_indices[0]] = vertex_normals[triangles[i].vertex_indices[0]] + v2;
		vertex_normals[triangles[i].vertex_indices[1]] = vertex_normals[triangles[i].vertex_indices[1]] + v2;
		vertex_normals[triangles[i].vertex_indices[2]] = vertex_normals[triangles[i].vertex_indices[2]] + v2;
	}

	for(size_t i = 0; i < vertex_normals.size(); i++)
		vertex_normals[i].normalize();
}

void mesh::generate_triangle_normals(void)
{
	if(triangles.size() == 0)
		return;

	triangle_normals.clear();
	triangle_normals.resize(triangles.size());

	for(size_t i = 0; i < triangles.size(); i++)
	{
		vertex_3 v0 = vertices[triangles[i].vertex_indices[1]] - vertices[triangles[i].vertex_indices[0]];
		vertex_3 v1 = vertices[triangles[i].vertex_indices[2]] - vertices[triangles[i].vertex_indices[0]];
		triangle_normals[i] = v0.cross(v1);
		triangle_normals[i].normalize();
	}
}

void mesh::generate_vertex_and_triangle_normals(void)
{
	generate_vertex_normals();
	generate_triangle_normals();
}

void mesh::regenerate_vertex_and_triangle_normals_if_exists(void)
{
	if(triangle_normals.size() > 0)
		generate_triangle_normals();

	if(vertex_normals.size() > 0)
		generate_vertex_normals();
}

template<typename T> void mesh::eliminate_vector_duplicates(vector<T> &v)
{
	if(0 == v.size())
		return;

	set<T> s(v.begin(), v.end()); // Eliminate duplicates (and sort them)
	vector<T> vtemp(s.begin(), s.end()); // Stuff things back into a temp vector
	v.swap(vtemp); // Assign temp vector contents to destination vector
}

bool mesh::merge_vertex_pair(const size_t keeper, const size_t goner)
{
	if(keeper >= vertices.size() || goner >= vertices.size())
		return false;

	if(keeper == goner)
		return true;

	// Merge vertex to triangle data.

	// Add goner's vertex to triangle data to keeper's triangle to vertex data,
	// and replace goner's index with keeper's index in all relevant triangles' index data.
	for(size_t i = 0; i < vertex_to_triangle_indices[goner].size(); i++)
	{
		size_t tri_index = vertex_to_triangle_indices[goner][i];

		vertex_to_triangle_indices[keeper].push_back(tri_index);

		for(size_t j = 0; j < 3; j++)
			if(goner == triangles[tri_index].vertex_indices[j])
				triangles[tri_index].vertex_indices[j] = keeper;
	}

	// Finalize keeper's vertex to triangle data.
	eliminate_vector_duplicates(vertex_to_triangle_indices[keeper]);

	// Clear out goner's vertex to triangle data for good.
	vertex_to_triangle_indices[goner].clear();


	// Merge vertex to vertex data.

	// Add goner's vertex to vertex data to keeper's vertex to vertex data,
	// and replace goner's index with keeper's index in all relevant vertices' vertex to vertex data.
	for(size_t i = 0; i < vertex_to_vertex_indices[goner].size(); i++)
	{
		size_t vert_index = vertex_to_vertex_indices[goner][i];

		vertex_to_vertex_indices[keeper].push_back(vert_index);

		for(size_t j = 0; j < vertex_to_vertex_indices[vert_index].size(); j++)
		{
			// Could probably break after this, but whatever...
			if(goner == vertex_to_vertex_indices[vert_index][j])
				vertex_to_vertex_indices[vert_index][j] = keeper;
		}

		eliminate_vector_duplicates(vertex_to_vertex_indices[vert_index]);
	}

	// Finalize keeper's vertex to vertex data.
	eliminate_vector_duplicates(vertex_to_vertex_indices[keeper]);

	// Clear out goner's vertex to vertex data for good.
	vertex_to_vertex_indices[goner].clear();

	// Note: At this point, vertices[goner] is now a rogue vertex.
	// We will skip erasing it from the vertices vector because that would mean a whole lot more work
	// (we'd have to reindex every vertex after it in the vector, etc.).
	// 
	// If the mesh is saved to STL, then the rogue vertex will automatically be skipped, and life is good.
	//
	// If the mesh is saved to POV-Ray mesh2, then the rogue vertex will be included in the vertex
	// list, but it will simply not be referenced in the triangle list -- this is a bit inoptimal
	// in terms of the file size (it will add a few dozen unneeded bytes to the file size).

	return true;
}
