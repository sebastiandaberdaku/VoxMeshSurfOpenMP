add_subdirectory(SeedFill3D)
add_subdirectory(PLY)
add_subdirectory(MC)

set(Geometry_source_files 
	HierarchicalQueue.h
	point3D.h
	voxel.h
	voxel_offset.cpp
	voxel_offset.h
	HierarchicalQueue.cpp
	point3D.cpp
	voxel.cpp
)

add_library(Geometry ${Geometry_source_files})
target_link_libraries(Geometry PLY SeedFill3D MC)
