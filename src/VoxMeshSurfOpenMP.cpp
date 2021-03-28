/*
 * TriangVoxSurfOpenMP.cpp
 *
 *  Created on: 9/jan/2018
 *      Author: Sebastian Daberdaku
 */
#include "Molecule/Molecule.h"
#include "CommandLineParser/CommandLineParser.h"
#include "utils/disclaimer.h"
#include "utils/makeDirectory.h"
#include <iostream>
#include <numeric>
#include <stdlib.h>
#include <string>
#include <utility>
#include "exceptions/ParsingPDBException.h"
using namespace std;

int main (int argc, char* argv[]) {
//---------------------------variables and parameters------------------------------------
	float solventRadius; // solvent-probe sphere radius

	float resolution; // resolution^3 = #voxels/Å^3
	string inname_molecule; // input filename
	string outname_molecule; // output filename
	string inname_radii; // input atom radii

	bool help = false; // if true, print the help message
	bool version = false; // if true, print the program version
	bool license = false; // if true, print the license information
	int  output_type;
	bool normalize_pose;
	bool hydrogen;
	bool hetatm;
	int smoothing_iterations;
	float lambda;
	float mu;
	bool vertex_normals;

	int surface_type;

    auto const t_start = std::chrono::high_resolution_clock::now();

	if (!CommandLineParser::parseCommandLine(argc, argv,
			solventRadius, resolution, surface_type,
			inname_molecule, outname_molecule, inname_radii,
			hydrogen, hetatm, normalize_pose, output_type,
			smoothing_iterations, lambda, mu, vertex_normals,
			help, version, license)) {
		return EXIT_FAILURE;
	}
	if (help || version || license) {
		return EXIT_SUCCESS;
	}
//-----------------------------print config-------------------------------------------

	PROGRAM_INFO
	/* summary of the parsed parameters */
	cout << "The specification is: \n";
	cout << "input filename: " << inname_molecule << "\n";
	cout << "output filename: " << outname_molecule << "\n";
	cout << "solvent-probe radius:\t" << solventRadius << "Å, \n";
	string sType;
	if (surface_type == 1)
		sType = "vdW";
	else if (surface_type == 2)
		sType = "SAS";
	else
		sType = "SES";

	cout << "surface type:\t" << sType << "\n";
	cout << "resolution:\t" << pow((double) resolution, 3.0) << " voxels per Å³, \n";
	if (normalize_pose)
		cout << "pose normalization: yes\n";
	else
		cout << "pose normalization: no\n";
	if (hetatm)
		cout << "include HETATM records: yes\n";
	else
		cout << "include HETATM records: no\n";
	if (hydrogen)
		cout << "include hydrogen atoms: yes\n";
	else
		cout << "include hydrogen atoms: no\n";
	cout << "atomic radii: " << inname_radii << "\n";
	cout << "**************************************************\n";

//-----------------------------computation--------------------------------------------
	try {
		Molecule m(surface_type, solventRadius, resolution, hydrogen, hetatm, normalize_pose,
				output_type, smoothing_iterations, lambda, mu, vertex_normals,
				inname_molecule, outname_molecule, inname_radii);
//-----------------------------conclusion --------------------------------------------
		cout << "**************************************************\n";
		cout << "Total calculation time:\t" << elapsedTime(t_start, std::chrono::high_resolution_clock::now()) << "\n";
		cout << "**************************************************\n";

	} catch (ParsingPDBException const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (fstream::failure const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (out_of_range const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (invalid_argument const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (logic_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	} catch (runtime_error const & e) {
		cerr << "error: " << e.what() << "\n";
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
