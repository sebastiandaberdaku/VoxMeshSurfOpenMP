/*
 * MolecularSurface.cpp
 *
 *  Created on: Dec 23, 2013
 *      Author: daberdaku
 */
/**
 * Implementation of the MolecularSurface class. This class
 * defines the molecular surface object, with all the methods
 * needed for its calculation.
 */
#include "../exceptions/ParsingPDBException.h"
#include "../Geometry/EDT/EDT.h"
#include "../Geometry/HierarchicalQueue.h"
#include "../Geometry/SeedFill3D/rapid3DSurfaceExtract.h"
#include "../Geometry/Sphere/DrawSphere.h"
#include "../utils/disclaimer.h"
#include "../utils/makeDirectory.h"
#include "../utils/numerical_limits.h"
#include "MolecularSurface.h"
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>

/**
 * Constructor of the class. Initializes the voxel grid and other data structures.
 *
 * \param atoms				List containing the atoms of the molecule.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * 							false: use original atomic radius values.
 * \param length			Length of the voxel grid.
 * \param width				Width of the voxel grid.
 * \param height			Height of the voxel grid.
 * \param translation		Translation vector to the grid's coordinate system
 * 							(already scaled).
 *
 */

MolecularSurface::MolecularSurface(int surface_type, int output_type, std::vector<atom> const & atoms,
		float solventRadius, float resolution, uint16_t length, uint16_t width, uint16_t height, point3D translation,
		bool normalize_pose, double const RI[3][3], point3D const & centroid) :
		atoms(&atoms), solventRadius(solventRadius), resolution(resolution),
		plength(length), pwidth(width), pheight(height), ptran(translation),
		surface(NULL), cpk_model(NULL), surface_type(surface_type), output_type(output_type),
		normalize_pose(normalize_pose), centroid(centroid){
	if (atoms.empty())
		throw invalid_argument("MolecularSurface::MolecularSurface() - The molecule has no atoms!");
	if (normalize_pose)
		for (int i = 0; i < 3; ++i)
			for (int j = 0; j < 3; ++j)
				this->RI[i][j] = RI[i][j];

} /* MolecularSurface() */

/**
 * Destructor of the class.
 */
MolecularSurface::~MolecularSurface() {
	if (surface != NULL)
		delete surface;
	if (cpk_model != NULL)
		delete cpk_model;
} /* ~MolecularSurface() */

/** Calculates a bounding box containing the molecule. This method
 * calculates the maximal and minimal coordinates reached by the
 * molecule by checking all the coordinates of the atoms composing it.
 *
 * \param atomsInModel		List containing all the model's atoms.
 * \param probeRadius		Desired probe radius of the solvent rolling sphere.
 * \param resolution		resolution^3 = number of voxels per Å^3
 * \param length			(return value) the length of the voxel grid.
 * \param width				(return value) the width of the voxel grid.
 * \param height			(return value) the height of the voxel grid.
 * \param translation		(return value) the translation vector to the
 * 							grid's coordinate system (already scaled).
 * \param max_atm_radii		(return value) the maximum atomic radii in
 * 							atomsInModel
 */
void MolecularSurface::boundingBox(std::vector<atom> const & atomsInModel,
		float solventRadius, float resolution, uint16_t & length,
		uint16_t & width, uint16_t & height, point3D & translation) {
	float max_atm_radius = min_float;
	if (atomsInModel.empty())
		throw ParsingPDBException("No atoms in PDB model.",
				"MolecularSurface::boundingBox", "Probably due to corrupt PDB file.");
	point3D minp(max_float, max_float, max_float);
	point3D maxp(min_float, min_float, min_float);
	for (auto const & atm : atomsInModel) {
		if (atm.x < minp.x)
			minp.x = atm.x;
		if (atm.y < minp.y)
			minp.y = atm.y;
		if (atm.z < minp.z)
			minp.z = atm.z;
		if (atm.x > maxp.x)
			maxp.x = atm.x;
		if (atm.y > maxp.y)
			maxp.y = atm.y;
		if (atm.z > maxp.z)
			maxp.z = atm.z;
		if (max_atm_radius < atm.radius)
			max_atm_radius = atm.radius;
	}
	float d = max_atm_radius + 2 * solventRadius + 1;

	maxp += point3D(d, d, d);
	minp -= point3D(d, d, d);

	/* transformation values */
	translation = -minp;

	/* bounding box dimensions */
	double boxLength = ceil(resolution * (maxp.x - minp.x)) + 2;
	double boxWidth  = ceil(resolution * (maxp.y - minp.y)) + 2;
	double boxHeight = ceil(resolution * (maxp.z - minp.z)) + 2;

	if ((boxLength <= UINT16_MAX)
			&& (boxWidth <= UINT16_MAX)
			&& (boxHeight <= UINT16_MAX)) {
		length = static_cast<uint16_t>(boxLength);
		width  = static_cast<uint16_t>(boxWidth);
		height = static_cast<uint16_t>(boxHeight);
		cout << "Grid size: " << length << " X " << width << " X " << height << endl;
	} else {
		std::stringstream ss;
		ss << "MolecularSurface::boundingBox() - ";
		ss << "The bounding box's dimensions exceed the maximum value of " << UINT16_MAX << ". ";
		ss << "Try setting a lower \"resolution\" value.";
		throw std::invalid_argument(ss.str());
	}
} /* boundingBox() */

/** Fills the voxels in the grid occupied by the molecule (protein).
 * This method implements a space-filling algorithm which is the
 * preliminary step for our grid-based macro-molecular surface
 * generation.
 *
 * In chemistry, a space-filling model, also known as a calotte model,
 * is a type of three-dimensional molecular model where the atoms are
 * represented by spheres whose radii are proportional to the radii of
 * the atoms and whose center-to-center distances are proportional to
 * the distances between the atomic nuclei, all in the same scale.
 * Atoms of different chemical elements are usually represented by
 * spheres of different colors.
 *
 * Calotte models are distinguished from other 3D representations,
 * such as the ball-and-stick and skeletal models, by the use of "full
 * size" balls for the atoms. They are useful for visualizing the
 * effective shape and relative dimensions of the molecule, in particular
 * the region of space occupied by it. On the other hand, calotte models
 * do not show explicitly the chemical bonds between the atoms, nor the
 * structure of the molecule beyond the first layer of atoms.
 *
 * Space-filling models are also called CPK models after the chemists
 * Robert Corey, Linus Pauling and Walter Koltun, who pioneered their use.
 */
void MolecularSurface::buildCPKModel() {
	if (cpk_model != NULL)
		delete cpk_model;
	cpk_model = new bool[plength * pwidth * pheight]();
//	unordered_map<int, vector<voxel_offset>> atom_sphere_offsets;

	float s = 0;
	if (surface_type > 1)
		s = solventRadius;
	/* For every atom in our list, calculate the voxels it occupies. */
#pragma omp parallel for schedule(runtime)
	for (size_t i = 0; i < atoms->size(); ++i) {
		/* Translate and discretize the coordinates */
		int cx = static_cast<int>(round(((*atoms)[i].x + ptran.x) * resolution));
		int cy = static_cast<int>(round(((*atoms)[i].y + ptran.y) * resolution));
		int cz = static_cast<int>(round(((*atoms)[i].z + ptran.z) * resolution));
		int r  = static_cast<int>(round(((*atoms)[i].radius + s ) * resolution));
		DrawBall(cpk_model, pwidth, pheight, cx, cy, cz, r);
//		DrawBallf(cpk_model, pwidth, pheight,
//				((*atoms)[i].x + ptran.x) * resolution,
//				((*atoms)[i].y + ptran.y) * resolution,
//				((*atoms)[i].z + ptran.z) * resolution,
//				((*atoms)[i].radius + s ) * resolution);

	}
} /* buildCPKModel() */


/**
 * Build the boundary of the solid created with the space-filling algorithm.
 * A voxel belongs to the boundary if it has at least one neighbor which is not
 * occupied by any atom.
 */
void MolecularSurface::buildSurface() {
	if (surface != NULL)
		delete surface;

	surface = new bool[plength * pwidth * pheight]();

	float s = 0;
	if (surface_type > 1)
		s = solventRadius;
	/*
	 * the model could have disjoint parts,
	 * all atom centers must be used as seeds
	 */
#pragma omp parallel for schedule(runtime)
	for (size_t i = 0; i < atoms->size(); ++i) {
		/* Translate and discretize the coordinates */
		int cx = static_cast<int>(round(((*atoms)[i].x + ptran.x) * resolution));
		int cy = static_cast<int>(round(((*atoms)[i].y + ptran.y) * resolution));
		int cz = static_cast<int>(round(((*atoms)[i].z + ptran.z) * resolution));
		int c  = static_cast<int>(round(((*atoms)[i].radius + s ) * resolution)) + 1;
		int d  = 2 * c + 1;
		bool * cpkm = new bool[d * d * d];
		bool * temp = new bool[d * d * d];
		bool * surf = new bool[d * d * d]();
		for (int ii = cx - c; ii <= cx + c; ++ii) {
			for (int jj = cy - c; jj <= cy + c; ++jj) {
				for (int kk = cz - c; kk <= cz + c; ++kk) {
					cpkm[((ii - cx + c) * d + jj - cy + c) * d + kk - cz + c] = cpk_model[(ii * pwidth + jj) * pheight + kk];
					temp[((ii - cx + c) * d + jj - cy + c) * d + kk - cz + c] = cpk_model[(ii * pwidth + jj) * pheight + kk];
				}
			}
		}
		rapid3DSurfaceExtract::surfaceExtract3D(cpkm, temp, surf, d, d, d, voxel(c, c, c));
		for (int ii = cx - c; ii <= cx + c; ++ii) {
			for (int jj = cy - c; jj <= cy + c; ++jj) {
				for (int kk = cz - c; kk <= cz + c; ++kk) {
					if (surf[((ii - cx + c) * d + jj - cy + c) * d + kk - cz + c])
						surface[(ii * pwidth + jj) * pheight + kk] = true;
				}
			}
		}
		delete[] temp;
		delete[] cpkm;
		delete[] surf;
	}
} /* buildSurface() */


void MolecularSurface::computeSAS() {
	cout << "Calculating Solvent Accessible CPK Model." << endl;
	buildCPKModel();
	if (output_type < 5) {
	    cout << "Building Solvent Accessible surface." << endl;
		double start_time = omp_get_wtime();
		buildSurface();
		cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	}
}

void MolecularSurface::computevdW() {
	cout << "Calculating van der Waals CPK Model." << endl;
	buildCPKModel();
	if (output_type < 5) {
	    cout << "Building van der Waals surface." << endl;
		double start_time = omp_get_wtime();
		buildSurface();
		cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	}

}

/** Calculates a partial Euclidean Distance Map of the voxelized
 * representation of the molecular surface. Starting from the Solvent-
 * Accessible Surface of the molecule, it is possible to obtain the
 * molecular surface (or Solvent-Excluded Surface) applying the EDT
 * (Euclidean Distance Transform). The SES is then extracted from the
 * partial Euclidean Distance Map, as it's voxels have a particular
 * distance value, i.e. greater or equal to the solvent/probe sphere
 * radius. The calculation starts from the boundary (surface) voxels,
 * and procedes towards the inside of the molecule, one shell at a time.
 * The initial shell is obviously the surface of the molecule. The
 * second shell is composed by voxels internal to the molecule, that
 * are neighbors to the first shell's voxels, and don't belong to the
 * first shell. The third shell will be composed by neighbors of the
 * second shell's voxels, internal to the molecule, and that don't
 * belong to the first or second shell, and so on. The distance map
 * values are propagated from the molecular surface towards the inside
 * of the molecule one shell at a time.
 */
void MolecularSurface::computeSES() {
	cout << "Calculating Solvent Accessible CPK Model." << endl;
	double start_time = omp_get_wtime();
	buildCPKModel();
	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	cout << "Building Solvent Accessible surface." << endl;
	start_time = omp_get_wtime();
	buildSurface();
	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;

	cout << "Calculating the Euclidean Distance Transform." << endl;

	start_time = omp_get_wtime();
	float min_sdist = powf(solventRadius * resolution, 2.0);

	uint16_t numberOfQueues = static_cast<uint16_t>(ceil(pow(solventRadius * resolution - 0.5, 2.0)));
	if (numberOfQueues < 2)
		numberOfQueues = 2;

	bool* temp_cpk_model = new bool[plength * pwidth * pheight];
#pragma omp parallel for schedule(runtime)
	for (size_t ii = 0; ii < plength * pwidth * pheight; ++ii)
		temp_cpk_model[ii] = cpk_model[ii];

	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;

	/* For every atom in our list, calculate the voxels it occupies. */
	start_time = omp_get_wtime();
#pragma omp parallel
{
		vector<uint16_t> dMap(plength * pwidth * pheight, UINT16_MAX);
		#pragma omp	for schedule(runtime)
		for (size_t ii = 0; ii < (*atoms).size(); ++ii) {
			computeEDT(temp_cpk_model, surface, numberOfQueues, min_sdist, (*atoms)[ii],
					ptran, resolution, solventRadius, pwidth, pheight, dMap, cpk_model);
		}
}
	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;

//#pragma omp parallel for schedule(runtime)
//	for (uint16_t ii = 0; ii < plength; ++ii)
//		for (uint16_t jj = 0; jj < pwidth; ++jj)
//			for (uint16_t kk = 0; kk < pheight; ++kk)
//			if (cpk_model[(ii * pwidth + jj) * pheight + kk]){
//				size_t nb = 0;
//				if (cpk_model[((ii+1) * pwidth + jj) * pheight + kk]) ++nb;
//				if (cpk_model[((ii-1) * pwidth + jj) * pheight + kk]) ++nb;
//				if (cpk_model[(ii * pwidth + jj+1) * pheight + kk]) ++nb;
//				if (cpk_model[(ii * pwidth + jj-1) * pheight + kk]) ++nb;
//				if (cpk_model[(ii * pwidth + jj) * pheight + kk+1]) ++nb;
//				if (cpk_model[(ii * pwidth + jj) * pheight + kk-1]) ++nb;
//				if (nb <= 1)
//					cpk_model[(ii * pwidth + jj) * pheight + kk] = false;
//			}

	if (output_type < 5) {
		cout << "Calculating Solvent Excluded Surface." << endl;
		start_time = omp_get_wtime();
		buildSurface();
		cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	}
	delete[] temp_cpk_model;
} /*  */

/** Prints the 3D voxelized representation to file using the PCD
 * (Point Cloud Data) file format.
 *
 * Each PCD file contains a header that identifies and declares
 * certain properties of the point cloud data stored in the file.
 * The header of a PCD must be encoded in ASCII.
 * Storing point cloud data in both a simple ascii form with each
 * point on a line, space or tab separated, without any other
 * characters on it, as well as in a binary dump format, allows
 * us to have the best of both worlds: simplicity and speed,
 * depending on the underlying application. The ascii format
 * allows users to open up point cloud files and plot them using
 * standard software tools like gnuplot or manipulate them using
 * tools like sed, awk, etc.
 *
 * For a detailed description of the PCD (Point Cloud Data) file
 * format specification see:
 * http://pointclouds.org/documentation/tutorials/pcd_file_format.php
 *
 * \param filename	Name of the output file. The '.pcd' extension
 * 					is added automatically.
 *
 * \throws ofstream::failure
 */
void MolecularSurface::outputPCDModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	string sType;
	if (surface_type == 1)
		sType = "vdW";
	else if (surface_type == 2)
		sType = "SAS";
	else
		sType = "SES";
	file_stream.open("./output/" + filename + "_" + sType + ".pcd");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(5);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (cpk_model[(i * pwidth + j) * pheight + k])
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	/* File header */
	file_stream << "# .PCD v.7 - Point Cloud Data file format\n"
			<< "# created with " << PROGRAM_NAME << "\n"
			<< "# version " << PROGRAM_VERSION << "\n"
			<< "VERSION .7\n" << "FIELDS x y z\n" << "SIZE 4 4 4\n"
			<< "TYPE F F F\n" << "COUNT 1 1 1\n" << "WIDTH "
			<< surfaceVoxels.size() << "\n" << "HEIGHT 1\n"
			<< "VIEWPOINT 0 0 0 1 0 0 0\n" << "POINTS " << surfaceVoxels.size()
			<< "\n" << "DATA ascii";

	while (!surfaceVoxels.empty()) {
		point3D currentVoxel(surfaceVoxels.front());
		currentVoxel /= resolution;
		currentVoxel -= ptran;
		if (normalize_pose) {
			currentVoxel = currentVoxel.transform(RI);
			currentVoxel += centroid;
		}
		file_stream << "\n" << currentVoxel.x
				<< " " << currentVoxel.y
				<< " " << currentVoxel.z;
		surfaceVoxels.pop();
	}
	file_stream.close();
} /* outputSurfacePCDModel() */



void MolecularSurface::outputOpenDXModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	string sType;
	if (surface_type == 1)
		sType = "vdW";
	else if (surface_type == 2)
		sType = "SAS";
	else
		sType = "SES";
	file_stream.open("./output/" + filename + "_" + sType + ".dx");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	size_t gridSize = plength * pwidth * pheight;

	file_stream  << std::setprecision(5);
	/* File header */
	file_stream << "# created with " << PROGRAM_NAME << "\n";
	file_stream << "# version " << PROGRAM_VERSION << "\n";
	file_stream << "object 1 class gridpositions counts " << plength << " "	<< pwidth << " " << pheight << "\n";
	file_stream << "origin " << -ptran.x << " " << -ptran.y << " " << -ptran.z << "\n";
	file_stream << "delta " << 1/resolution << " " << 0.0f << " " << 0.0f << "\n";
	file_stream << "delta " << 0.0f << " " << 1/resolution << " " << 0.0f << "\n";
	file_stream << "delta " << 0.0f << " " << 0.0f << " " << 1/resolution << "\n";
	file_stream << "object 2 class gridpositions counts " << plength << " "	<< pwidth << " " << pheight << "\n";
	file_stream << "object 3 class array type double rank 0 items "<< gridSize << " data follows\n";
	/* data */
	size_t c = 1;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface[(i * pwidth + j) * pheight + k])
					file_stream << 1.0f;
				else
					file_stream << 0.0f;
				if (c % 3 == 0 || c == gridSize) {
					file_stream << "\n";
				}
				else {
					file_stream << " ";
				}
				++c;
			}
		}
	}
	file_stream << "attribute \"dep\" string \"positions\"\n";
	file_stream	<< "object \"regular positions regular connections\" class field\n";
	file_stream << "component \"positions\" value 1\n";
	file_stream << "component \"connections\" value 2\n";
	file_stream << "component \"data\" value 3\n";
	file_stream.close();
} /* outputSurfaceOpenDXModel() */

void MolecularSurface::outputVTKStructuredPointsModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	string sType;
	if (surface_type == 1)
		sType = "vdW";
	else if (surface_type == 2)
		sType = "SAS";
	else
		sType = "SES";
	file_stream.open("./output/" + filename + "_" + sType + ".vtk");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}

	file_stream  << std::setprecision(5);
	/* File header */
	file_stream << "# vtk DataFile Version 3.0\n";
	file_stream << filename << ".vtk created with " << PROGRAM_NAME << " - version " << PROGRAM_VERSION << ".\n";
	file_stream << "ASCII\n";
	file_stream << "DATASET STRUCTURED_POINTS\n";
	file_stream << "DIMENSIONS " << plength << " " << pwidth << " " << pheight << "\n";
	file_stream << "ORIGIN " << -ptran.x << " " << -ptran.y << " " << -ptran.z << "\n";
	file_stream << "SPACING " << 1/resolution<< " " << 1/resolution<< " " << 1/resolution<< "\n";
	file_stream << "POINT_DATA " << plength * pwidth * pheight << "\n";
	file_stream << "SCALARS " << filename << " int 1\n";
	file_stream << "LOOKUP_TABLE default";

	/* data */
	for (int k = 0; k < pheight; ++k) {
		for (int j = 0; j < pwidth; ++j) {
			for (int i = 0; i < plength; ++i) {
				file_stream << "\n";
				if (surface[(i * pwidth + j) * pheight + k])
					file_stream << "1";
				else
					file_stream << "0";
			}
		}
	}

	file_stream.close();
} /* outputVTKStructuredPointsModel() */

void MolecularSurface::outputVTKPolyDataModel(std::string const & filename) {
	std::ofstream file_stream;
	makeDirectory("./output");
	string sType;
	if (surface_type == 1)
		sType = "vdW";
	else if (surface_type == 2)
		sType = "SAS";
	else
		sType = "SES";
	file_stream.open("./output/" + filename + "_" + sType + ".vtk");
	if (!file_stream.is_open()) {
		throw std::ofstream::failure("Error opening output file.");
	}
	file_stream  << std::setprecision(5);
	std::queue<voxel> surfaceVoxels;
	for (int i = 0; i < plength; ++i) {
		for (int j = 0; j < pwidth; ++j) {
			for (int k = 0; k < pheight; ++k) {
				if (surface[(i * pwidth + j) * pheight + k])
					surfaceVoxels.push(voxel(i, j, k));
			}
		}
	}
	/* File header */
	file_stream << "# vtk DataFile Version 3.0\n";
	file_stream << filename << ".vtk created with " << PROGRAM_NAME << " - version " << PROGRAM_VERSION <<".\n";
	file_stream << "ASCII\n";
	file_stream << "DATASET POLYDATA\n";
	file_stream << "POINTS " << surfaceVoxels.size() << " float";

	/* data */
	while (!surfaceVoxels.empty()) {
		point3D currentVoxel(surfaceVoxels.front());
		currentVoxel /= resolution;
		currentVoxel -= ptran;
		if (normalize_pose) {
			currentVoxel = currentVoxel.transform(RI);
			currentVoxel += centroid;
		}
		file_stream << "\n" << currentVoxel.x
				<< " " << currentVoxel.y
				<< " " << currentVoxel.z;
		surfaceVoxels.pop();
	}

	file_stream.close();
} /* outputSurfaceVTKModel() */


void MolecularSurface::runMarchingCubes(int smoothing_iterations, float lambda, float mu, bool vertex_normals) {
	double start_time = omp_get_wtime();
//	marchingCubes(cpk_model, plength, pwidth, pheight, faces, vertices);
	marchingCubesParallel(cpk_model, plength, pwidth, pheight, faces, vertices);

	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	cout << "Total of number of faces: " << faces.size() << endl;
	cout << "Total of number of vertices: " << vertices.size() << endl;

	if (smoothing_iterations > 0) {
		cout << "Laplacian smoothing... " << endl;
		start_time = omp_get_wtime();
		LaplacianSmoothing2(smoothing_iterations, lambda, mu, faces, vertices, resolution);
		cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	}
	cout << "Scaling and translating mesh... " << endl;
	start_time = omp_get_wtime();
#pragma omp parallel for schedule(runtime)
	for (size_t ii = 0; ii < vertices.size(); ++ii) {
		vertices[ii] /= resolution;
		vertices[ii] -= ptran;
		if (normalize_pose) {
			vertices[ii] = vertices[ii].transform(RI);
			vertices[ii] += centroid;
		}
	}
	if (vertex_normals) {
		cout << "Computing normals... " << endl;
		start_time = omp_get_wtime();
		computeNormals(faces, vertices, normals);
		cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
	}
	cout << "elapsed time: " << omp_get_wtime() - start_time << endl;
}

void MolecularSurface::outputPLY(string const & filename, bool binary) {
	using namespace tinyply;
	makeDirectory("./output");
	std::filebuf fb;
	if (binary)
		fb.open("./output/" + filename + "-binary.ply", std::ios::out | std::ios::binary);
	else
		fb.open("./output/" + filename + "-ascii.ply", std::ios::out);
	std::ostream outstream_ply(&fb);
	outstream_ply << std::setprecision(10);
	if (outstream_ply.fail())
		throw std::runtime_error("failed to open " + filename);

	PlyFile output_mesh;
	output_mesh.add_properties_to_element("vertex", { "x", "y", "z" },
		Type::FLOAT64, vertices.size(), reinterpret_cast<uint8_t*>(vertices.data()), Type::INVALID, 0);
	if (!normals.empty())
		output_mesh.add_properties_to_element("vertex", { "nx", "ny", "nz" },
				Type::FLOAT64, normals.size(), reinterpret_cast<uint8_t*>(normals.data()), Type::INVALID, 0);
	output_mesh.add_properties_to_element("face", { "vertex_indices" },
		Type::UINT32, faces.size(), reinterpret_cast<uint8_t*>(faces.data()), Type::UINT8, 3);

	stringstream ss;
	ss << "PLY file created with " << PROGRAM_NAME << " - version " << PROGRAM_VERSION;
    output_mesh.get_comments().push_back(ss.str());
	output_mesh.write(outstream_ply, binary);
}
