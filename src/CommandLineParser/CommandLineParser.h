/*
 * CommandLineParser.h
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

#ifndef COMMANDLINEPARSER_H_
#define COMMANDLINEPARSER_H_

#include <string>

using namespace std;

/**
 * This class provides a simple interface for accessing the command line arguments.
 */
class CommandLineParser {
public:
	static bool parseCommandLine(int const & argc, char const * const argv[],
			float & solventRadius, float & resolution, int & surface_type,
			string & inname, string & outname, string & inname_radii,
			bool & hydrogen, bool & hetatm, bool & normalize_pose,
			int & output_type, int & smoothing_iterations, float & lambda,
			float & mu, bool & vertex_normals, bool & help, bool & version,
			bool & license);
};

#endif /* COMMANDLINEPARSER_H_ */
