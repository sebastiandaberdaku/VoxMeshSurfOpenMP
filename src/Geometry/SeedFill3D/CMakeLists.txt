set(SeedFill3D_source_files 
	node.h
	range.h
	sliceNode.cpp
	node.cpp
	range.cpp
	rapid3DSurfaceExtract.h
	sliceNode.h
)

add_library(SeedFill3D ${SeedFill3D_source_files})
