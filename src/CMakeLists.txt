cmake_minimum_required (VERSION 3.1 FATAL_ERROR)
project (VoxMeshSurfOpenMP)

set(CMAKE_CXX_COMPILER g++)

include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
	message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()

find_package(OpenMP)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O2")

find_package(Boost COMPONENTS regex system program_options thread REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

add_definitions(-DNDEBUG -D_GLIBCXX_PARALLEL) #-D NO_OUTPUT_TEST

add_subdirectory(PDB)
add_subdirectory(Molecule)
add_subdirectory(MolecularSurface)
add_subdirectory(Geometry)
add_subdirectory(exceptions)
add_subdirectory(CommandLineParser)

set(VoxMeshSurfOpenMP_source_files VoxMeshSurfOpenMP.cpp)
add_executable(VoxMeshSurfOpenMP ${VoxMeshSurfOpenMP_source_files})
target_link_libraries(VoxMeshSurfOpenMP Molecule exceptions CommandLineParser ${Boost_LIBRARIES}) 
