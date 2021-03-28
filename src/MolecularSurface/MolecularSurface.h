/*
 * MolecularSurface.h
 *
 *  Created on: Nov 11, 2013
 *      Author: daberdaku
 */
/**
 * This header defines the MolecularSurface object, representing
 * the surface of a molecule, with the related methods for its
 * calculation.
 */
#ifndef MOLECULARSURFACE_H_
#define MOLECULARSURFACE_H_

#include "../Geometry/point3D.h"
#include "../PDB/atom.h"

#include <algorithm>
#include <functional>

#include <list>
#include <queue>
#include <vector>
#include <omp.h>
#include <unordered_map>
#include "../Geometry/MC/LaplacianSmoothing.h"
#include "../Geometry/PLY/tinyply.h"

using namespace std;

class MolecularSurface {
	friend class Molecule;
public:
	bool *surface, *cpk_model;
	float resolution;

	MolecularSurface(int surface_type, int output_type, vector<atom> const & atoms, float solventRadius,
			float resolution, uint16_t length, uint16_t width, uint16_t height, point3D translation,
			bool normalize_pose, double const RI[3][3], point3D const & centroid);

	virtual ~MolecularSurface();

	static void boundingBox(vector<atom> const & atomsInModel,
			float probeRadius, float resolution, uint16_t & length,
			uint16_t & width, uint16_t & height, point3D & translation);

	void computeSES();
	void computeSAS();
	void computevdW();

	void outputPCDModel(std::string const & filename);
	void outputOpenDXModel(string const & filename);
	void outputVTKStructuredPointsModel(string const & filename);
	void outputVTKPolyDataModel(string const & filename);
	void outputPLY(string const & filename, bool binary = true);

	void runMarchingCubes(int smoothing_iterations, float lambda, float mu, bool vertex_normals);

private:
	void buildCPKModel();
	void buildSurface();

	static inline void reduce_distance_maps(vector<uint16_t> * output,
			vector<uint16_t> const * input) {
		transform(output->begin(), output->end(), input->begin(), output->begin(),
				[](uint16_t a, uint16_t b) {return min(a, b);} );
	}

#pragma omp declare reduction(dmap_reduction : vector<uint16_t> * : reduce_distance_maps(omp_out, omp_in)) initializer(omp_priv(omp_orig))

	static inline void reduce_usets(unordered_set<triangle> & output,
			unordered_set<triangle> const & input) {
		output.insert(input.begin(), input.end());
	}

#pragma omp declare reduction(uset_reduction : unordered_set<triangle> : reduce_usets(omp_out, omp_in)) initializer(omp_priv(omp_orig))

	point3D ptran;
	vector<atom> const * atoms;
	float solventRadius;
	uint16_t pheight, pwidth, plength;
	int surface_type, output_type;

	vector<faceinfo> faces;
	vector<point3D> vertices;
	vector<point3D> normals;

	bool normalize_pose;
	point3D centroid;
	double RI[3][3];
};
#endif /* MOLECULARSURFACE_H_ */
