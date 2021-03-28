/*
 * voxel.h
 *
 *  Created on: 29/ago/2014
 *      Author: sebastian
 */

#ifndef VOXEL_H_
#define VOXEL_H_

#include "voxel_offset.h"
#include <inttypes.h>
#include <ostream>
/**
 * Struct defining a voxel data type.
 * A single voxel is univocally identified by it's coordinates.
 */
#if defined TIGHT_PACKING
#pragma pack(push)  /* push current alignment to stack */
#pragma pack(1)     /* set alignment to 1 byte boundary */
#endif
typedef struct voxel {
	int32_t ix, iy, iz;
	voxel();
	voxel(int32_t i, int32_t j, int32_t k);
	voxel(voxel const & vox);
//	voxel & operator=(voxel const & vox);
	inline int32_t sDistance_to(voxel const & x) const {
		int32_t dx, dy, dz;
		dx = ix - x.ix;
		dy = iy - x.iy;
		dz = iz - x.iz;
		return (dx*dx + dy*dy + dz*dz);
	}
	friend inline std::ostream & operator<<(std::ostream& os, voxel const & v) {
		os << "[" << v.ix << ", " << v.iy << ", " << v.iz << "]";
		return os;
	};

	inline voxel & operator+=(voxel_offset const & rhs) {
		this->ix += rhs.i;
		this->iy += rhs.j;
		this->iz += rhs.k;
		return *this;
	};
	inline const voxel operator+(voxel_offset const & rhs) const {
		return voxel(*this) += rhs;
	};
	inline voxel & operator-=(voxel_offset const & rhs) {
		this->ix -= rhs.i;
		this->iy -= rhs.j;
		this->iz -= rhs.k;
		return *this;
	};
	inline const voxel operator-(voxel_offset const & rhs) const {
		return voxel(*this) -= rhs;
	};
	inline const voxel_offset operator-(voxel const & rhs) const {
		int32_t i = this->ix - rhs.ix;
		int32_t j = this->iy - rhs.iy;
		int32_t k = this->iz - rhs.iz;
		return voxel_offset(i, j, k);
	};
	/**
	 * Copy assignment operator
	 */
	inline voxel & operator=(voxel const & vox) {
		if (this != &vox) {
			this->ix = vox.ix;
			this->iy = vox.iy;
			this->iz = vox.iz;
		}
		return *this;
	}
	inline bool operator==(voxel const & rhs) const {
		return (this->ix == rhs.ix && this->iy == rhs.iy && this->iz == rhs.iz);
	};

} voxel;
#if defined TIGHT_PACKING
#pragma pack(pop)   /* restore original alignment from stack */
#endif

#include <boost/functional/hash.hpp>
namespace std {

template<>
struct hash<voxel> {
	size_t operator()(const voxel& t) const {
		// Start with a hash value of 29    .
		size_t seed = 29;
		// Modify 'seed' by XORing in order to get the same hash
		// if the three vertices are the same but shuffled among them.
		boost::hash_combine(seed, hash<int32_t>()(t.ix));
		boost::hash_combine(seed, hash<int32_t>()(t.iy));
		boost::hash_combine(seed, hash<int32_t>()(t.iz));
		// Return the result.
		return seed;
	}
};
}

#endif /* VOXEL_H_ */
