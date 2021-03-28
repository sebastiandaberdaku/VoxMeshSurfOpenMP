#include "primitives.h"
#include "mesh.h"

#include <string>
using std::string;


int main(int argc, char **argv)
{
	if(argc != 2)
	{
		cout << "Example usage: tsmoothstl filename.stl" << endl;
		return 1;
	}

	mesh mesh;

	if(false == mesh.load_from_binary_stereo_lithography_file(argv[1]))
	{
		cout << "Error: Could not properly read file " << argv[1] << endl;
		return 2;
	}


	mesh.fix_cracks();


	// Set mesh size via different properties

//	mesh.set_max_extent(3);
//	mesh.set_area(3);
//	mesh.set_volume(3);


	// Smooth mesh via different schemes

//	float lambda = 0.5f;
//	float mu = -0.53f;
//	mesh.taubin_smooth(lambda, mu, 20);
//	mesh.taubin_smooth_fujiwara(lambda, mu, 20);
//	mesh.taubin_smooth_cn(lambda, mu, 20);


	string out_file_name = argv[1];
	out_file_name = "processed_" + out_file_name;

	if(false == mesh.save_to_binary_stereo_lithography_file(out_file_name.c_str()))
	{
		cout << "Error: Could not properly write file " << out_file_name << endl;
		return 2;
	}


	// Save to POV-Ray mesh2 format
//	mesh.save_to_povray_mesh2_file("filename.inc");


	return 0;
}

