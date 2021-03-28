set(Molecule_source_files
	Molecule.h
	Molecule.cpp
)

add_library(Molecule ${Molecule_source_files})
target_link_libraries(Molecule MolecularSurface)
