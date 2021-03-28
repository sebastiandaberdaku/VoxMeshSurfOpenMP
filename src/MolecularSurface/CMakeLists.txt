set(MolecularSurface_source_files
	MolecularSurface.h
	MolecularSurface.cpp
)

add_library(MolecularSurface ${MolecularSurface_source_files})
target_link_libraries(MolecularSurface Geometry PDB)
