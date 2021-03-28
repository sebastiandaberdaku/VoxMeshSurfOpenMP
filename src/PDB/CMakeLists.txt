set(PDB_source_files
	atom.h
	atom.cpp
	PDBModel.h
	PDBModel.cpp
	charmm22.h
)

add_library(PDB ${PDB_source_files})
