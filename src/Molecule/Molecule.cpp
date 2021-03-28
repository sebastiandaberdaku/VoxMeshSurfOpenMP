/*
 * Molecule.cpp
 *
 *  Created on: 07/set/2014
 *      Author: sebastian
 */

#include "../utils/disclaimer.h"
#include "../utils/doubleCompare.h"
#include "../utils/makeDirectory.h"
#include "Molecule.h"
#include <unordered_map>

using namespace std;

Molecule::Molecule(int surface_type, float solventRadius,float resolution,
		bool hydrogen, bool hetatm, bool normalize_pose, int output_type,
		int smoothing_iterations, float lambda, float mu, bool vertex_normals,
		string const & inname, string const & outname,
		string const & inname_radii) :
	centroid(point3D(0, 0, 0)), inname_molecule(inname), resolution(resolution) {

	cout << "Loading PDB file: "<< inname <<"\n";
	pdbModel = new PDBModel(inname, inname_radii, hydrogen, hetatm);
	surface = NULL;

	if (normalize_pose) {
		poseNormalization(pdbModel->atomsInModel, atoms);
		invert(R, RI);
	}
	else
		atoms = pdbModel->atomsInModel;

	sort(begin(atoms), end(atoms), [](atom const & a, atom const & b) {return a.x < b.x;});


    double t_start = omp_get_wtime();
    cout << "Initializing parameters for " << inname << ".\n";
	cout << "Number of atoms in model: " << atoms.size() << "\n";

	cout << "Number of parallel threads: " << omp_get_max_threads() << "\n";

	omp_sched_t schedule_type;
	int chunk_size;

	omp_get_schedule(&schedule_type, &chunk_size);

	string policy[] = {"static", "dynamic", "guided", "auto"};
	cout << "OpenMP loop scheduling policy: " << policy[schedule_type-1] << " chunk size: " << chunk_size << "\n";

	MolecularSurface::boundingBox(atoms, solventRadius, resolution, length,
			width, height, translation);

	if (surface != NULL)
		delete surface;
	surface = new MolecularSurface(surface_type, output_type, atoms, solventRadius,
			resolution, length, width, height, translation,
			normalize_pose, RI, centroid);
	switch (surface_type) {
	case 1:
		cout << "Calculating van der Waals surface.\n";
		surface->computevdW();
		break;
	case 2:
		cout << "Calculating Solvent Accessible surface.\n";
		surface->computeSAS();
		break;
	default:
		cout << "Calculating Solvent Excluded surface.\n";
		surface->computeSES();
	}

	switch (output_type) {
	case 1:
		cout << "Output surface PCD model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputPCDModel(outname);
#endif
		break;
	case 2:
		cout << "Output surface OpenDX model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputOpenDXModel(outname);
#endif
		break;
	case 3:
		cout << "Output surface VTK (Structured Points) model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputVTKStructuredPointsModel(outname);
#endif
		break;
	case 4:
		cout << "Output surface VTK (PolyData) model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputVTKPolyDataModel(outname);
#endif
		break;
	case 5:
		cout << "Running marching cubes algorithm.\n";
		surface->runMarchingCubes(smoothing_iterations, lambda, mu, vertex_normals);
		cout << "Output surface mesh (ASCII PLY) model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputPLY(outname, false);
#endif
		break;
	default :
		cout << "Running marching cubes algorithm.\n";
		surface->runMarchingCubes(smoothing_iterations, lambda, mu, vertex_normals);
		cout << "Output surface mesh (binary PLY) model.\n";
#ifndef NO_OUTPUT_TEST
		surface->outputPLY(outname, true);
#endif
	}

	cout << "Total computation time (in seconds):\t" << omp_get_wtime() - t_start << "\n";
}

Molecule::~Molecule() {
	if (surface != NULL)
		delete surface;
	if (pdbModel != NULL)
		delete pdbModel;
}

