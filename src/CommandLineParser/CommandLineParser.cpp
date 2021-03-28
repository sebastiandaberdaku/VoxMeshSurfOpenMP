/*
 * CommandLineParser.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#include "../utils/disclaimer.h"
#include "CommandLineParser.h"
#include "CustomValidators.h"

/**
 * This method parses the input command line options.
 * @param argc		input argc reference, i.e. the first argument of the main method
 * @param argv		input argv reference, i.e. the second argument of the main method
 *
 * @return 			true if all arguments are parsed correctly, false otherwise.
 */
bool CommandLineParser::parseCommandLine(int const & argc,
		char const * const argv[], float & solventRadius, float & resolution,
		int & surface_type, string & inname, string & outname,
		string & inname_radii, bool & hydrogen, bool & hetatm,
		bool & normalize_pose, int & output_type, int & smoothing_iterations,
		float & lambda, float & mu, bool & vertex_normals, bool & help,
		bool & version, bool & license) {

	//------------------------------- command line options --------------------------------------
	stringstream ss;
	ss << PROGRAM_NAME << " - version " << PROGRAM_VERSION << ". Usage";
	options_description description(ss.str());
	description.add_options()("help,h", "Display this help message.")
	/*
	 * This token is used to specify the input filename. It is mandatory.
	 * The pqr command line token has no option name. The command line tokens which have
	 * no option name are called "positional options" in the program_options library.
	 */
	("pdb", value<input_PDB_filename>()->required(), "Input PDB or XYZR file.")
	/*
	 * This token is used to specify the output filename.
	 */
	("output,o", value<filename>(), "Output filename. If not specified, the input filename (with no extension) will be used.")
	/*
	 * This token is used to specify the input atom radii filename.
	 */
	("atom_radii", value<filename>(), "File containing the radius information of each atom. If not specified, the default CHARMM22 radius values will be used.")
	/*
	 * This token is used to specify the output type.
	 */
	("output_type,t", value<out_type>()->default_value(6), "Output type (default is 6):\n1-Point Cloud Data file (*.pcd);\n2-OpenDX scalar data format (*.dx);\n3-Visualization Toolkit Structured Points (*.vtk);\n4-Visualization Toolkit PolyData (*.vtk);\n5-Surface mesh ASCII model (*.ply);\n6-Surface mesh binary model (*.ply).")
	("surface_type,s", value<surf_type>()->default_value(3), "Surface type (default is 3):\n1-van der Waals;\n2-Solvent Accessible Surface;\n3-Solvent Excluded Surface.")
	/*
	 * The -p (--probe_radius) flag is used to specify the probe radius. Default is 1.4. Optional.
	 */
	("probe_radius,p", value<probe_radius>()->default_value(1.4), "Solvent-probe radius (in Å), positive floating point number (default is 1.4).")
	/*
	 * The -r (--resolution) flag is used to specify the resolution of the voxel grid. Default is 4. Optional.
	 */
	("resolution,r", value<resolution_param>()->default_value(2.0), "Resolution factor, positive floating point (default is 2.0). This value's cube determines the number of voxels per Å³.")
	("Laplacian_smoothing,L", value<num_iterations>()->default_value(5), "Number of Laplacian smoothing iterations on the output mesh (default is 1). If set to 0, no smoothing will be carried out on the surface mesh.")
	("lambda", value<shrink_param>()->default_value(0.90), "lambda parameter in the Laplacian smoothing step (positive floating point, default is 0.90).")
	("mu", value<inflate_param>()->default_value(-0.92), "mu parameter in the Laplacian smoothing step (negative floating point, default is -0.92).")
	("no_vertex_normals", "Skip the vertex normals estimation step. The output mesh will not have smooth shading.")
	/*
	 * The --hetatm flag is used to include HETATM records in the surface computation.
	 */
	("hetatm", "Include HETATM records in the surface computation (only for *.pdb input files, has no effect for *.xyzr files).")
	/*
	 * The --hetatm flag is used to include HETATM records in the surface computation.
	 */
	("hydrogen", "Include hydrogen atoms in the surface computation (only for *.pdb input files, has no effect for *.xyzr files).\nNOTE: X-ray crystallography cannot resolve hydrogen atoms in most protein crystals, so in most PDB files, hydrogen atoms are absent. Sometimes hydrogens are added by modeling. Hydrogens are always present in PDB files resulting from NMR analysis, and are usually present in theoretical models.")
	/*
	 * The --normalize_pose flag is used to run pose normalization on the molecule
	 * information.
	 */
	("no_pose_normalization", "By default, the three principal axes of the molecule are aligned with the three Cartesian axes in order to reduce the size of the voxel grid containing the molecule. If this option is set, no pose normalization will take place.")
	/*
	 * The (--license) flag is used to view the program's license
	 * information.
	 */
	("license", "View license information.")
	/*
	 * The -v (--version) flag is used to view the program's version
	 * information.
	 */
	("version,v", "Display the version number");
	/*
	 * The input filename must be declared as positional.
	 */
	positional_options_description p;
	p.add("pdb", 1);

	variables_map vm;
	try {
		//--------------------------------parsing command line options------------------------------
		/*
		 * And it is finally specified when parsing the command line.
		 */
		store(command_line_parser(argc, argv).options(description).positional(p).run(), vm);

		if (vm.count("help")) {
			cout << description;
			help = true;
			return true;
		}
		if (vm.count("version")) {
			cout << "Program: " << argv[0] << ", version: " << PROGRAM_VERSION << "\n";
			version = true;
			return true;
		}
		if (vm.count("license")) {
			DISCLAIMER
			license = true;
			return true;
		}
		if (vm.count("no_pose_normalization"))
			normalize_pose = false;
		else
			normalize_pose = true;
		if (vm.count("hetatm"))
			hetatm = true;
		else
			hetatm = false;
		if (vm.count("hydrogen"))
			hydrogen = true;
		else
			hydrogen = false;
		if (vm.count("no_vertex_normals"))
			vertex_normals = false;
		else
			vertex_normals = true;

		/*
		 * notify throws exceptions so we call it after the above checks
		 */
		notify(vm);

		/* initializing variables */
		inname = vm["pdb"].as<input_PDB_filename>().filename;
		if (vm.count("output")) {
			outname = vm["output"].as<filename>().fname;
		} else {
			int lastindex = inname.find_last_of(".");
			outname = inname.substr(0, lastindex);
			lastindex = outname.find_last_of("/");
			outname = outname.substr(lastindex + 1);
		}
		if (vm.count("atom_radii")) {
			inname_radii = vm["atom_radii"].as<filename>().fname;
		} else {
			inname_radii = "";
		}

		smoothing_iterations = vm["Laplacian_smoothing"].as<num_iterations>().iter;
		lambda = vm["lambda"].as<shrink_param>().lambda;
		mu = vm["mu"].as<inflate_param>().mu;
		solventRadius = vm["probe_radius"].as<probe_radius>().p;
		surface_type = vm["surface_type"].as<surf_type>().st;


		resolution = vm["resolution"].as<resolution_param>().r;
		output_type = vm["output_type"].as<out_type>().st;

	} catch (error const & e) {
		cerr << "error: " << e.what() << "\n";
		cerr << description << "\n";
		return false;
	}
	return true;
};


