/*
 * EDT.h
 *
 *  Created on: 17/mar/2018
 *      Author: Sebastian Daberdaku
 */

#ifndef EDT_EDT_H_
#define EDT_EDT_H_

#include "../../PDB/atom.h"
#include "../HierarchicalQueue.h"
#include "../voxel_offset.h"
#include "../voxel.h"

inline void computeEDT(bool const * cpk_model, bool const * sa_surface, uint16_t numberOfQueues,
		float min_sdist, atom const & atm, point3D const & ptran, float resolution,
		float probeRadius, uint16_t pwidth, uint16_t pheight, vector<uint16_t> & distanceMap, bool * new_cpk_model){

	HierarchicalQueue* HQ1 = new HierarchicalQueue(numberOfQueues);
	HierarchicalQueue* HQ2 = new HierarchicalQueue(numberOfQueues);

	int cx = static_cast<int>(round((atm.x + ptran.x) * resolution));
	int cy = static_cast<int>(round((atm.y + ptran.y) * resolution));
	int cz = static_cast<int>(round((atm.z + ptran.z) * resolution));
	uint16_t d = static_cast<int>(round((atm.radius +  probeRadius) * resolution)) + 3;
	uint16_t e = static_cast<int>(round(probeRadius * resolution)) + 1;
	uint16_t L = 2 * d + 2 * e + 1 ;
	vector<voxel> nearestSurfaceVoxel(L * L * L);

	int ii, jj, kk;
	/* for all voxels */
	for (int i = cx - d; i <= cx + d; ++i) {
		ii = i - cx + d + e;
		for (int j = cy - d; j <= cy + d; ++j) {
			jj = j - cy + d + e;
			for (int k = cz - d; k <= cz + d; ++k) {
				kk = k - cz + d + e;
				/* create the first shell by pushing the surface voxels in the list */
				if (sa_surface[(i * pwidth + j) * pheight + k]) {
					HQ1->push(voxel(i, j, k), 0);
					/* a surface voxel has zero distance from itself */
					nearestSurfaceVoxel[(ii * L + jj) * L + kk] = voxel(i, j, k);
					distanceMap[(i * pwidth + j) * pheight + k] = 0;
					new_cpk_model[(i * pwidth + j) * pheight + k] = false;
				}
			}
		}
	}
	while (!HQ1->empty()) {
		/*current voxel*/
		voxel cv = HQ1->front();
		HQ1->pop();
		ii = cv.ix + d + e - cx;
		jj = cv.iy + d + e - cy;
		kk = cv.iz + d + e - cz;
		voxel nearestSurfVox = nearestSurfaceVoxel[(ii * L + jj) * L + kk];
		uint16_t squaredDistance = distanceMap[(cv.ix * pwidth + cv.iy) * pheight + cv.iz];
		voxel nb; // neighbour
		bool isEnd = true;
		for (int i = 0; i < 26; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);
				ii = nb.ix + d + e - cx;
				jj = nb.iy + d + e - cy;
				kk = nb.iz + d + e - cz;
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
						&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(ii * L + jj) * L + kk] = nearestSurfVox;
					HQ1->push(nb, newSDistance);
					isEnd = false;

					if (newSDistance < min_sdist)
						new_cpk_model[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = false;
				}
			}
		}
		if (isEnd && squaredDistance >= 24) {
			HQ2->push(cv, squaredDistance);
		}
	}

	delete HQ1;
	HQ1 = NULL;

	while (!HQ2->empty()) {
		/*current voxel*/
		voxel cv = HQ2->front();
		HQ2->pop();

		ii = cv.ix - cx + d + e;
		jj = cv.iy - cy + d + e;
		kk = cv.iz - cz + d + e;
		voxel nearestSurfVox = nearestSurfaceVoxel[(ii * L + jj) * L + kk];// neighbor
		voxel nb;
		for (int i = 0; i < 124; ++i) {
			nb = cv + voxel_offset::neighbours[i];
			if (cpk_model[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]) {
				uint16_t newSDistance = nb.sDistance_to(nearestSurfVox);

				ii = nb.ix - cx + d + e;
				jj = nb.iy - cy + d + e;
				kk = nb.iz - cz + d + e;
				if (newSDistance < distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz]
					&& newSDistance < numberOfQueues) {
					distanceMap[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = newSDistance;
					nearestSurfaceVoxel[(ii * L + jj) * L + kk] = nearestSurfVox;
					HQ2->push(nb, newSDistance);

					if (newSDistance < min_sdist)
						new_cpk_model[(nb.ix * pwidth + nb.iy) * pheight + nb.iz] = false;
				}
			}
		}
	}
	delete HQ2;
	HQ2 = NULL;
}

#endif /* EDT_EDT_H_ */
