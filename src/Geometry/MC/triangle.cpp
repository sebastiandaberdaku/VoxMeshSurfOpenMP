/*
 * triangle.cpp
 *
 *  Created on: Nov 15, 2018
 *      Author: sebastian
 */



#include "triangle.h"

triangle::triangle() :
		p1(0, 0, 0), p2(0, 0, 0), p3(0, 0, 0) { };
triangle::triangle(double x1, double y1, double z1,
			double x2, double y2, double z2,
			double x3, double y3, double z3) :
		p1(x1, y1, z1), p2(x2, y2, z2), p3(x3, y3, z3) { };
triangle::triangle(point3D const & p1, point3D const & p2, point3D const & p3) :
		p1(p1), p2(p2), p3(p3) { };
triangle::triangle(triangle const & t) :
		p1(t.p1), p2(t.p2), p3(t.p3) { };
