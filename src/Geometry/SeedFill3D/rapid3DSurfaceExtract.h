/*
 * rapid3DSurfaceExtract.h
 *
 *  Created on: Feb 5, 2015
 *      Author: sebastian
 */
/**
 * Here we implement a method that extracts the molecular surface from
 * a volumetric model, using the 3D seed filling-algorithm proposed in:
 * Wei-Wei Yu, Fei He, Ping Xi, A rapid 3D seed-filling algorithm based on scan slice,
 * Computers & Graphics, Volume 34, Issue 4, August 2010, Pages 449-459, ISSN 0097-8493,
 * http://dx.doi.org/10.1016/j.cag.2010.05.005.
 * http://www.sciencedirect.com/science/article/pii/S0097849310000725
 */
#ifndef SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_
#define SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_

#include "node.h"
#include "sliceNode.h"
#include <stack>
#include "../voxel_offset.h"
#include <array>

namespace rapid3DSurfaceExtract {

/**
 * Clears a whole range in temp, from xl to xr, on line y of slice z.
 * Also checks if any of the voxels in the range belongs to the surface.
 * @param image		target image
 * @param temp		current 3D image
 * @param output	output voxelGrid
 * @param zl		range lower bound
 * @param zr		range upper bound
 * @param y			range line
 * @param x			image slice
 */
static inline void clearRange(bool const * imageData, bool * temp, bool * output, uint16_t length, uint16_t width, uint16_t height,
		int zl, int zr, int x, int y) {
	for (int k = zl; k <= zr; k++) {
		if ((x > 0) && (x < length - 1) && (y > 0) && (y < width - 1)) {
			if (!imageData[(x * width + y + 1) * height + k]
			   ||!imageData[(x * width + y - 1) * height + k]
			   ||!imageData[((x + 1) * width + y) * height + k]
			   ||!imageData[((x - 1) * width + y) * height + k]
//			   ||!imageData[((x + 1) * width + y + 1) * height + k]
//			   ||!imageData[((x + 1) * width + y - 1) * height + k]
//			   ||!imageData[((x - 1) * width + y + 1) * height + k]
//			   ||!imageData[((x - 1) * width + y - 1) * height + k]
							)
				output[(x * width + y) * height + k] = true;
		}
		temp[(x * width + y) * height + k] = false;
	}
	if (zl > 0)
		output[(x * width + y) * height + zl] = true;
	if (zr < length - 1)
		output[(x * width + y) * height + zr] = true;
};

/**
 * Returns the filled and valid range discovered scanning the target image from seed
 * @param temp		current 3D image
 * @param seed		the seed
 * OUTPUT:
 * @param zl 		the left range delimiter
 * @param zr 		the right range delimiter
 */
static inline void obtainFilledValidRange(bool const * temp, uint16_t length, uint16_t width, uint16_t height,
		voxel const & seed, int & zl, int & zr) {
	zr = seed.iz;
	while ((zr < height) && temp[(seed.ix * width + seed.iy) * height + zr]) {
		++zr;
	}
	--zr;
	zl = seed.iz;
	while ((zl > -1) && temp[(seed.ix * width + seed.iy) * height + zl]) {
		--zl;
	}
	++zl;
};

/**
 * Checks if the current range [xpl, xpr] contains any valid seeds. If a seed is found
 * the corresponding range is extracted, filled and inserted in dList.
 * @param imageData	the current voxelGrid
 * @param temp		temporary voxelGrid
 * @param output	output voxelGrid
 * @param zpl		lower bound of the range to be checked
 * @param zpr		upper bound of the range to be checked
 * @param y			line of the current range
 * @param x			slice of the current range
 * OUTPUT:
 * @param dList		list of the filled ranges in the current slice
 * @param zl		lower bound of the newly found range
 * @param zr		upper bound of the newly found range
 * @return			true if a valid seed is found, false otherwise
 */
static inline bool checkAndClearRange(bool const * imageData, bool * temp, bool * output,
		uint16_t length, uint16_t width, uint16_t height,
		int zpl, int zpr, int x, int y, IdList & dList, int & zl, int & zr) {
	for (int k = zpl; k <= zpr; ++k) {
		if (temp[(x * width + y) * height + k]) { // if a valid seed is found
			// extract the unfilled range,and write the filled seeds ID into dList
			// search the valid and unfilled range ([zl,zr]) from [zpl, zpr]
			obtainFilledValidRange(temp, length, width, height, voxel(x, y, k), zl, zr);
			//extract and fill the range [zl,zr], mark the seeds and write the ID of them into dList
			clearRange(imageData, temp, output, length, width, height, zl, zr, x, y);
			dList.push_back(range(zl, zr, y));
			return true;
		}
	}
	return false;
};

/**
 * This method flood-fills a single slice of the input voxelGrid, starting from voxel seed.
 * The slice is identified by the z coordinate of seed. The ID of each seed is extracted from
 * dList.
 * @param imageData	target voxelGrid
 * @param temp		temporary voxelGrid
 * @param output	output voxelGrid
 * @param seed		starting voxel for the flood-filling procedure
 * OUTPUT:
 * @param seedList	list of the filled ranges in the current slice
 */
static inline void surfaceExtract2D(bool const * imageData, bool * temp, bool * output,
		uint16_t length, uint16_t width, uint16_t height, voxel const & seed, IdList & seedList) {
	// empty the seed list
	seedList.clear();
	// initialize the empty 2D stack
	stack<node> stack2D;
	// obtain the unfilled and valid range ([zl,zr]) by using seed
	int zl, zr; // range [zl, zr]
	obtainFilledValidRange(temp, length, width, height, seed, zl, zr);
	// extract and fill the range [zl,zr],
	clearRange(imageData, temp, output, length, width, height,  zl, zr, seed.ix, seed.iy);
	// mark the seeds and write the ID of them into dList
	seedList.push_back(range(zl, zr, seed.iy));
	// create two nodes with range ([zl, zr]), push them into 2D stack, along two opposite directions
	if (seed.iy > 0)
		stack2D.push(node(zl, zr, seed.iy, -1));
	if (seed.iy < width - 1)
		stack2D.push(node(zl, zr, seed.iy, 1));
	// obtain the length of this range
	while (!stack2D.empty()) { // if the 2D stack is not empty
		//pop a node and then set the search range ([xpl,xpr]) onto the next scan line
		node c(stack2D.top()); // current node
		stack2D.pop();
		int zpl = c.zl, zpr = c.zr, y = c.y + c.direction_y;

		while (checkAndClearRange(imageData, temp, output, length, width, height, zpl, zpr, seed.ix, y, seedList, zl, zr)) { //if the range [xpl, xpr] is valid
			if (zr < zpr - 1) { //it maybe exist other valid and unfilled ranges on the same scan line
				for (int k = zr + 2; k <= zpr; ++k) { // search the leftmost seed of the next valid range, which locates at the same scan line
					if (temp[(seed.ix * width + y) * height + k]) { /* a valid seed exists */
						stack2D.push(node(zr + 2, zpr, y - c.direction_y, c.direction_y));
						break;
					}
				}
			}
			//execute the necessary rollback operation
			if (zl < zpl - 1) {	//the rollback operation may occur on the left side of xpl
			// rollback to the previous scan line
				int new_y = y - c.direction_y;
				// search a valid and unfilled seed on the left side of xpl
				for (int k = zl; k <= zpl - 2; ++k) {
					if (temp[(seed.ix * width + new_y) * height + k]) { /* a valid seed exists */
						//the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [xl, xpl-2], the direction is opposite
						stack2D.push(node(zl, zpl - 2, y, -c.direction_y));
						break;
					}
				}
			}
			if (zr > zpr + 1) {	//the rollback operation may occur on the right side of xpr
			// rollback to the previous scan line
				int new_y = y - c.direction_y;
				// search a valid and unfilled seed on the right side of xpr
				for (int k = zpr + 2; k <= zr; ++k) {
					// if (a valid seed exists){//the search range is valid, push it into 2D stack
					if (temp[(seed.ix * width + new_y) * height + k]) { /* a valid seed exists */
						// the search range is valid, push it into 2D stack
						// push a new node into 2D stack, the range is [xpr+ 2, xr], the direction is opposite
						stack2D.push(node(zpr + 2, zr, y, -c.direction_y));
						break;
					}
				}
			}
			// continue the loop from this range [xl, xr], along the same direction as before
			// obtain the new search range by assigning the range [xl, xr] to [xpl, xpr]
			zpl = zl;
			zpr = zr;
			y += c.direction_y;
			if (y < 0 || y >= width)
				break;
		}
	}
};

/**
 * This method projects the seed list of the current slice onto slice x and returns true if
 * a valid seed is found in the projected area, false otherwise.
 * @param temp		the current voxelGrid
 * @param seed_list	the list of seeds in the current slice
 * @param x			the neighbor slice
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 * OUTPUT:
 * @param sy 		y coordinate of the first valid seed found
 * @param sz		z coordinate of the first valid seed found
 * @return			true if a new seed is found, false otherwise
 */
static inline bool searchNeighborSlice(bool const * temp, uint16_t length, uint16_t width, uint16_t height, IdList & seed_list,
		int x, int leap_var, int & sy, int & sz) {
	IdList::iterator id, end;
	for (id = seed_list.begin(), end = seed_list.end(); id != end; ++id) {
		for (int iz = id->zl; iz <= id->zr; iz += leap_var) {
			if (temp[(x * width + id->y) * height + iz]) {
				sz = iz;
				sy = id->y;
				/*
				 * Erase previous seeds, excluding the current one!
				 */
				seed_list.erase(seed_list.begin(), id);
				/*
				 * Resize the current seed, as all voxels from zl to the current iz are occupied
				 */
				id->zl = iz;
				return true;
			}
		}
	}
	seed_list.clear();
	return false;
};

/**
 * Here we implement a rapid 3D seed-filling algorithm, to extract the portion
 * of surface of the object in imageData which is connected to the given seed.
 * @param imageData input image
 * @param temp		temporary voxel grid, initially a copy of imageData
 * @param output	output voxelGrid
 * @param seed		starting voxel
 * @param leap_var	the leaping variable is used to leap over plenty of invalid
 * 					voxels in order to accelerate search process for the first
 * 					valid seed
 */
static inline void surfaceExtract3D(bool const * imageData, bool * temp, bool * output,
		uint16_t length, uint16_t width, uint16_t height, voxel const & seed, int leap_var = 1) {
	if (!imageData[(seed.ix * width + seed.iy) * height + seed.iz])
		return;
	//create empty 3D stack
	stack<sliceNode> stack3D;

	IdList currentSeeds;
	//extract and fill the initial unfilled region by using an initial seed (seed)
	surfaceExtract2D(imageData, temp, output, length, width, height, seed, currentSeeds);
	//create two nodes with this initial region(dList), push them into 3D stack, along two opposite directions
	if (seed.ix < length - 1)
		stack3D.push(sliceNode(currentSeeds, seed.ix, +1));
	if (seed.ix > 0)
		stack3D.push(sliceNode(currentSeeds, seed.ix, -1));
	while (!stack3D.empty()) {
		//pop a node from 3D stack
		sliceNode n(stack3D.top());
		stack3D.pop();
		//set the search region (PassList) by assigning the region (dList) to it
		currentSeeds = n.seedList; /*
							 * std::list::operator= c++11 Assigns new contents to the container,
							 * replacing its current contents, and modifying its size accordingly.
							 */
		int new_x = n.x + n.direction_x;
		//obtain a search region (Plist_mapped) by mapping the region (PassList) onto the next slice
		//and then check the validation of this mapped region (Plist_mapped)
		int sy, sz; // new seed coordinates
		while(searchNeighborSlice(temp, length, width, height, currentSeeds, new_x, leap_var, sy, sz)){ //if the region is valid
			//set a value to the leaping variable
			IdList newSeeds;
			//search a valid seed in the mapped region (Plist_mapped)
			//extract and fill this new unfilled region (newSeeds here represents the first new region)
			surfaceExtract2D(imageData, temp, output, length, width, height, voxel(new_x, sy, sz), newSeeds);
			//continue to search other new unfilled regions in the region (Plist_mapped)
			while (searchNeighborSlice(temp, length, width, height, currentSeeds, new_x, leap_var, sy, sz)){
				//extract and fill these new unfilled regions
				IdList otherSeeds;
				surfaceExtract2D(imageData, temp, output, length, width, height, voxel(new_x, sy, sz), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_x < length - 1)
					stack3D.push(sliceNode(otherSeeds, new_x, +1));
				if (new_x > 0)
					stack3D.push(sliceNode(otherSeeds, new_x, -1));
			}
			//execute the rollback operation
			//map the region (CurrentList) onto the previous slice
			currentSeeds = newSeeds;
			while (searchNeighborSlice(temp, length, width, height, newSeeds, new_x - n.direction_x, leap_var, sy, sz)) {
				IdList otherSeeds;
				//extract and fill these new unfilled regions
				surfaceExtract2D(imageData, temp, output, length, width, height,
						voxel(new_x - n.direction_x, sy, sz), otherSeeds);
				//create two nodes with each region, push them into 3D stack, along two opposite directions
				if (new_x - n.direction_x < length - 1)
					stack3D.push(sliceNode(otherSeeds, new_x - n.direction_x, +1));
				if (new_x - n.direction_x > 0)
					stack3D.push(sliceNode(otherSeeds, new_x - n.direction_x, -1));
			}
			//continue the loop from the region (CurrentList), along the same direction as before
			//obtain the new search region by assigning the region (CurrentList) to (PassList)
			new_x += n.direction_x;
			if (new_x < 0 || new_x == length)
				break;
		}
	}
};
}
#endif /* SEEDFILL3D_RAPID3DSURFACEEXTRACT_H_ */
