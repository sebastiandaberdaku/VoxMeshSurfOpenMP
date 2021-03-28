/*
 * triangle.h
 *
 *  Created on: Oct 25, 2018
 *      Author: sebastian
 */

#ifndef TRIANGLE_H_
#define TRIANGLE_H_

#include "../point3D.h"
#include <boost/functional/hash.hpp>

typedef struct triangle{
public:
	point3D p1, p2, p3;
	triangle();
	triangle(double x1, double y1, double z1,
			double x2, double y2, double z2,
			double x3, double y3, double z3);
	triangle(point3D const & p1, point3D const & p2, point3D const & p3);
	triangle(triangle const & t);
	/**
	 * Copy assignment operator
	 */
	inline triangle & operator=(triangle const & t) {
		if (this != &t) {
			this->p1 = t.p1;
			this->p2 = t.p2;
			this->p3 = t.p3;
		}
		return *this;
	};

	inline bool operator==(triangle const & rhs) const {
		return (this->p1 == rhs.p1 && this->p2 == rhs.p2 && this->p3 == rhs.p3)
				|| (this->p1 == rhs.p1 && this->p2 == rhs.p3 && this->p3 == rhs.p2)
				|| (this->p1 == rhs.p2 && this->p2 == rhs.p1 && this->p3 == rhs.p3)
				|| (this->p1 == rhs.p2 && this->p2 == rhs.p3 && this->p3 == rhs.p1)
				|| (this->p1 == rhs.p3 && this->p2 == rhs.p1 && this->p3 == rhs.p2)
				|| (this->p1 == rhs.p3 && this->p2 == rhs.p2 && this->p3 == rhs.p1);
	};

} triangle;

namespace std {

template<>
struct hash<triangle> {
	std::size_t operator()(const triangle& t) const {
		// Start with a hash value of 0    .
		std::size_t seed = 37;
		// Modify 'seed' by XORing in order to get the same hash
		// if the three vertices are the same but shuffled among them.
		seed ^= hash<point3D>()(t.p1);
		seed ^= hash<point3D>()(t.p2);
		seed ^= hash<point3D>()(t.p3);
		// Return the result.
		return seed;
	}
};
}

#endif /* TRIANGLE_H_ */
