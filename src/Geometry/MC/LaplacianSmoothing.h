/*
 * LaplacianSmoothing.h
 *
 *  Created on: Nov 7, 2018
 *      Author: sebastian
 */

#ifndef LAPLACIANSMOOTHING_H_
#define LAPLACIANSMOOTHING_H_

#include <unordered_set>
#include "marchingCubes.h"

using namespace std;

static inline void reduce_point3D_vector(vector<point3D> & output, vector<point3D> const & input) {
	for (size_t ii = 0; ii < output.size(); ++ii) {
		output[ii] += input[ii];
	}
}

#pragma omp declare reduction(point3D_vec_reduction : vector<point3D> : reduce_point3D_vector(omp_out, omp_in)) initializer(omp_priv(omp_orig))

static inline void computeNormals(vector<faceinfo> const & faces, vector<point3D> const & vertices,
		vector<point3D> & normals) {
	normals = vector<point3D>(vertices.size(), point3D(0, 0, 0));

#pragma omp parallel for schedule(runtime) reduction(point3D_vec_reduction:normals)
	for (size_t ii = 0; ii < faces.size(); ++ii) {
		point3D const & p1 = vertices[faces[ii].a];
		point3D const & p2 = vertices[faces[ii].b];
		point3D const & p3 = vertices[faces[ii].c];

		point3D n = (p1 - p2) * (p1 - p3);
		normals[faces[ii].a] += n;
		normals[faces[ii].b] += n;
		normals[faces[ii].c] += n;
	}
#pragma omp parallel for schedule(runtime)
	for (size_t ii = 0; ii < normals.size(); ++ii)
		normals[ii].normalize();
}


static inline void reduce_vertex_nb(vector<vector<uint32_t>> & output, vector<vector<uint32_t>> const & input) {
	size_t vertnumber = output[0].size();
	for (size_t ii = 0; ii < vertnumber; ++ii) {
		size_t outvert = output[0][ii];
		for (size_t j = 0; j < input[0][ii]; ++j) {
			bool newvert = true;
			for (size_t k = 0; k < outvert; ++k) {
				if (input[j+1][ii] == output[k+1][ii]) {
					newvert = false;
					break;
				}
			}
			if (newvert) {
				output[0][ii]++;
				output[output[0][ii]][ii] = input[j+1][ii];
			}
		}
	}

}

#pragma omp declare reduction(vertex_nb_reduction : vector<vector<uint32_t>> : reduce_vertex_nb(omp_out, omp_in)) initializer(omp_priv(omp_orig))

static inline void oneRingNeighbourhood(vector<faceinfo> const & faces, vector<point3D> const & vertices,
		vector<vector<uint32_t>> & vertex_nb) {
	size_t vertnumber = vertices.size();
	size_t facenumber = faces.size();
	const int numthreads = omp_get_max_threads();
	vector<vector<vector<uint32_t>>> vertex_nb_th(numthreads);

//#pragma omp parallel for schedule(runtime) reduction(vertex_nb_reduction:vertex_nb)
#pragma omp parallel
{
	const int ithread = omp_get_thread_num();
	vertex_nb_th[ithread] = vector<vector<uint32_t>>(20, vector<uint32_t>(vertnumber, 0));

#pragma omp for schedule(runtime)
	for (size_t i = 0; i < facenumber; i++) {
		// insert neighbours of vertex a
		bool newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].a]; ++j) {
			if (faces[i].b == vertex_nb_th[ithread][j+1][faces[i].a]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].a]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].a]][faces[i].a] = faces[i].b;
		}
		newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].a]; ++j) {
			if (faces[i].c == vertex_nb_th[ithread][j+1][faces[i].a]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].a]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].a]][faces[i].a] = faces[i].c;
		}
		// insert neighbours of vertex b
		newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].b]; ++j) {
			if (faces[i].a == vertex_nb_th[ithread][j+1][faces[i].b]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].b]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].b]][faces[i].b] = faces[i].a;
		}
		newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].b]; ++j) {
			if (faces[i].c == vertex_nb_th[ithread][j+1][faces[i].b]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].b]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].b]][faces[i].b] = faces[i].c;
		}
		// insert neighbours of vertex c
		newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].c]; ++j) {
			if (faces[i].a == vertex_nb_th[ithread][j+1][faces[i].c]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].c]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].c]][faces[i].c] = faces[i].a;
		}
		newvert = true;
		for (size_t j = 0; j < vertex_nb_th[ithread][0][faces[i].c]; ++j) {
			if (faces[i].b == vertex_nb_th[ithread][j+1][faces[i].c]) {
				newvert = false;
				break;
			}
		}
		if (newvert) {
			vertex_nb_th[ithread][0][faces[i].c]++;
			vertex_nb_th[ithread][vertex_nb_th[ithread][0][faces[i].c]][faces[i].c] = faces[i].b;
		}
	}
}

vertex_nb = vertex_nb_th[0];
#pragma omp parallel for schedule(runtime)
for (size_t ii = 0; ii < vertnumber; ++ii) {
	for (size_t ithread = 1; ithread < numthreads; ++ithread) {
		size_t outvert = vertex_nb[0][ii];
		for (size_t j = 0; j < vertex_nb_th[ithread][0][ii]; ++j) {
			bool newvert = true;
			for (size_t k = 0; k < outvert; ++k) {
				if (vertex_nb_th[ithread][j+1][ii] == vertex_nb[k+1][ii]) {
					newvert = false;
					break;
				}
			}
			if (newvert) {
				vertex_nb[0][ii]++;
				vertex_nb[vertex_nb[0][ii]][ii] = vertex_nb_th[ithread][j+1][ii];
			}
		}
	}
}

}



static inline void LaplacianSmoothing(int smoothing_iterations, float lambda, float mu,
		vector<faceinfo> & faces, vector<point3D> & vertices,
		float resolution) {
	size_t vertnumber = vertices.size();
	size_t facenumber = faces.size();

	/**
	 * vertex_nb is a 2D array that indicates the immediate neighbours of all vertices in the mesh.
	 * vertex_nb[0][i] contains the number of neighbours for vertex i, i.e. vertices[i]
	 * say vertex_nb[0][i] is N, then the indices of the neighbours of vertices[i] are
	 * in vertex_nb[1][i], vertex_nb[2][i], ..., vertex_nb[N][i]
	 */
	vector<vector<uint32_t>> vertex_nb;

	oneRingNeighbourhood(faces, vertices, vertex_nb);

	for (size_t k = 0; k < smoothing_iterations; k++) {
		vector<point3D> tps(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (uint32_t j = 0; j < vertex_nb[0][i]; ++j) {
				tps[i] += vertices[vertex_nb[j+1][i]];
			}
			tps[i] /= vertex_nb[0][i];
			tps[i] -= vertices[i];
			tps[i] *= lambda;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (uint32_t j = 0; j < vertex_nb[0][i]; ++j) {
				tps[i] += vertices[vertex_nb[j+1][i]];
			}
			tps[i] /= vertex_nb[0][i];
			tps[i] -= vertices[i];
			tps[i] *= mu;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
	}

}

static inline void LaplacianSmoothing2(int smoothing_iterations, float lambda, float mu,
		vector<faceinfo> & faces, vector<point3D> & vertices,
		float resolution) {
	size_t vertnumber = vertices.size();
	size_t facenumber = faces.size();

	/**
	 * vertex_nb is a 2D array that indicates the immediate neighbours of all vertices in the mesh.
	 * vertex_nb[0][i] contains the number of neighbours for vertex i, i.e. vertices[i]
	 * say vertex_nb[0][i] is N, then the indices of the neighbours of vertices[i] are
	 * in vertex_nb[1][i], vertex_nb[2][i], ..., vertex_nb[N][i]
	 */
	vector<vector<uint32_t>> vertex_nb;

	oneRingNeighbourhood(faces, vertices, vertex_nb);

	for (size_t k = 0; k < smoothing_iterations; k++) {
		vector<point3D> tps(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (uint32_t j = 0; j < vertex_nb[0][i]; ++j) {
				tps[i] += vertices[vertex_nb[j + 1][i]];
			}
			tps[i] /= vertex_nb[0][i];
			tps[i] -= vertices[i];
			tps[i] *= lambda;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
		tps = vector<point3D>(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			if (vertex_nb[0][i] > 5) {
				for (uint32_t j = 0; j < vertex_nb[0][i]; ++j) {
					tps[i] += vertices[vertex_nb[j + 1][i]];
				}
				tps[i] /= vertex_nb[0][i];
				tps[i] -= vertices[i];
				tps[i] *= mu;
			}
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
	}

}

static inline void LaplacianSmoothing3(int smoothing_iterations, float lambda, float mu,
		vector<faceinfo> & faces, vector<point3D> & vertices,
		float resolution) {
	size_t vertnumber = vertices.size();
	size_t facenumber = faces.size();

	/**
	 * vertex_nb is a 2D array that indicates the immediate neighbours of all vertices in the mesh.
	 * vertex_nb[0][i] contains the number of neighbours for vertex i, i.e. vertices[i]
	 * say vertex_nb[0][i] is N, then the indices of the neighbours of vertices[i] are
	 * in vertex_nb[1][i], vertex_nb[2][i], ..., vertex_nb[N][i]
	 */
	vector<vector<uint32_t>> vertex_nb;

	oneRingNeighbourhood(faces, vertices, vertex_nb);

	for (size_t k = 0; k < smoothing_iterations; k++) {
		vector<point3D> tps(vertnumber, point3D(0, 0, 0));
		vector<unordered_set<uint32_t>> twoRingNeighbourhood(vertnumber);

		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (uint32_t j = 1; j <= vertex_nb[0][i]; ++j) {
				for (uint32_t k = 1; k <= vertex_nb[0][vertex_nb[j][i]]; ++k) {
					twoRingNeighbourhood[i].insert(vertex_nb[k][vertex_nb[j][i]]);
				}
			}
			for (auto const &  j : twoRingNeighbourhood[i]) {
				tps[i] += vertices[j];
			}
			tps[i] /= twoRingNeighbourhood[i].size();
			tps[i] -= vertices[i];
			tps[i] *= lambda;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
		tps = vector<point3D>(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (auto const &  j : twoRingNeighbourhood[i]) {
				tps[i] += vertices[j];
			}
			tps[i] /= twoRingNeighbourhood[i].size();
			tps[i] -= vertices[i];
			tps[i] *= mu;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
	}

}

static inline void newSmoothing(int smoothing_iterations, float beta,
		vector<faceinfo> & faces, vector<point3D> & vertices,
		float resolution) {
	size_t vertnumber = vertices.size();
	size_t facenumber = faces.size();

	/**
	 * vertex_nb is a 2D array that indicates the immediate neighbours of all vertices in the mesh.
	 * vertex_nb[0][i] contains the number of neighbours for vertex i, i.e. vertices[i]
	 * say vertex_nb[0][i] is N, then the indices of the neighbours of vertices[i] are
	 * in vertex_nb[1][i], vertex_nb[2][i], ..., vertex_nb[N][i]
	 */
	vector<vector<uint32_t>> vertex_nb;

	oneRingNeighbourhood(faces, vertices, vertex_nb);

	for (size_t k = 0; k < smoothing_iterations; k++) {
		vector<point3D> tps(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			for (uint32_t j = 1; j <= vertex_nb[0][i]; ++j) {
				tps[i] += vertices[vertex_nb[j][i]];
			}
			tps[i] /= vertex_nb[0][i];
			tps[i] -= vertices[i];
			if (vertex_nb[0][i] > 4)
				tps[i] *= beta;
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
		tps = vector<point3D>(vertnumber, point3D(0, 0, 0));
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			if (vertex_nb[0][i] > 4) {
			for (uint32_t j = 1; j <= vertex_nb[0][i]; ++j) {
				tps[i] += vertices[vertex_nb[j][i]];
			}
			tps[i] /= vertex_nb[0][i];
			tps[i] -= vertices[i];
			tps[i] *= beta;
			}
		}
		#pragma omp parallel for schedule(runtime)
		for (size_t i = 0; i < vertnumber; i++) {
			vertices[i] += tps[i];
		}
	}

}

#endif /* LAPLACIANSMOOTHING_H_ */
