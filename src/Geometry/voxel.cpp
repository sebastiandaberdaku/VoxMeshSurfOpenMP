/*
 * voxel.cpp
 *
 *  Created on: Sep 6, 2014
 *      Author: sebastian
 */

/**
 * Implementation of the voxel object.
 */
#include "voxel.h"

/**
 * Default constructor
 */
voxel::voxel() :
		ix(0), iy(0), iz(0) { }
/**
 * Constructor
 */
voxel::voxel(int32_t i, int32_t j, int32_t k) :
		ix(i), iy(j), iz(k) { }
/**
 * Copy constructor.
 */
voxel::voxel(voxel const & vox) :
		ix(vox.ix), iy(vox.iy), iz(vox.iz) { }

