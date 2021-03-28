/*
 * mcubes.h
 *
 *  Created on: Oct 18, 2018
 *      Author: sebastian
 */

#ifndef MARCHINGCUBES_H_
#define MARCHINGCUBES_H_

#include "../../PDB/atom.h"
#include "../voxel.h"
#include "triangle.h"
#include <array>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <stdlib.h>
#include <unordered_map>
#include <vector>

using namespace std;

typedef struct faceinfo {
	uint32_t a, b, c;
} faceinfo;

static const int edgeTable[256]={
	0x0  , 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
	0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00,
	0x190, 0x99 , 0x393, 0x29a, 0x596, 0x49f, 0x795, 0x69c,
	0x99c, 0x895, 0xb9f, 0xa96, 0xd9a, 0xc93, 0xf99, 0xe90,
	0x230, 0x339, 0x33 , 0x13a, 0x636, 0x73f, 0x435, 0x53c,
	0xa3c, 0xb35, 0x83f, 0x936, 0xe3a, 0xf33, 0xc39, 0xd30,
	0x3a0, 0x2a9, 0x1a3, 0xaa , 0x7a6, 0x6af, 0x5a5, 0x4ac,
	0xbac, 0xaa5, 0x9af, 0x8a6, 0xfaa, 0xea3, 0xda9, 0xca0,
	0x460, 0x569, 0x663, 0x76a, 0x66 , 0x16f, 0x265, 0x36c,
	0xc6c, 0xd65, 0xe6f, 0xf66, 0x86a, 0x963, 0xa69, 0xb60,
	0x5f0, 0x4f9, 0x7f3, 0x6fa, 0x1f6, 0xff , 0x3f5, 0x2fc,
	0xdfc, 0xcf5, 0xfff, 0xef6, 0x9fa, 0x8f3, 0xbf9, 0xaf0,
	0x650, 0x759, 0x453, 0x55a, 0x256, 0x35f, 0x55 , 0x15c,
	0xe5c, 0xf55, 0xc5f, 0xd56, 0xa5a, 0xb53, 0x859, 0x950,
	0x7c0, 0x6c9, 0x5c3, 0x4ca, 0x3c6, 0x2cf, 0x1c5, 0xcc ,
	0xfcc, 0xec5, 0xdcf, 0xcc6, 0xbca, 0xac3, 0x9c9, 0x8c0,
	0x8c0, 0x9c9, 0xac3, 0xbca, 0xcc6, 0xdcf, 0xec5, 0xfcc,
	0xcc , 0x1c5, 0x2cf, 0x3c6, 0x4ca, 0x5c3, 0x6c9, 0x7c0,
	0x950, 0x859, 0xb53, 0xa5a, 0xd56, 0xc5f, 0xf55, 0xe5c,
	0x15c, 0x55 , 0x35f, 0x256, 0x55a, 0x453, 0x759, 0x650,
	0xaf0, 0xbf9, 0x8f3, 0x9fa, 0xef6, 0xfff, 0xcf5, 0xdfc,
	0x2fc, 0x3f5, 0xff , 0x1f6, 0x6fa, 0x7f3, 0x4f9, 0x5f0,
	0xb60, 0xa69, 0x963, 0x86a, 0xf66, 0xe6f, 0xd65, 0xc6c,
	0x36c, 0x265, 0x16f, 0x66 , 0x76a, 0x663, 0x569, 0x460,
	0xca0, 0xda9, 0xea3, 0xfaa, 0x8a6, 0x9af, 0xaa5, 0xbac,
	0x4ac, 0x5a5, 0x6af, 0x7a6, 0xaa , 0x1a3, 0x2a9, 0x3a0,
	0xd30, 0xc39, 0xf33, 0xe3a, 0x936, 0x83f, 0xb35, 0xa3c,
	0x53c, 0x435, 0x73f, 0x636, 0x13a, 0x33 , 0x339, 0x230,
	0xe90, 0xf99, 0xc93, 0xd9a, 0xa96, 0xb9f, 0x895, 0x99c,
	0x69c, 0x795, 0x49f, 0x596, 0x29a, 0x393, 0x99 , 0x190,
	0xf00, 0xe09, 0xd03, 0xc0a, 0xb06, 0xa0f, 0x905, 0x80c,
	0x70c, 0x605, 0x50f, 0x406, 0x30a, 0x203, 0x109, 0x0
};

static const int triTable[256][16] ={
	{	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	8,	3,	9,	8,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	1,	2,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	2,	10,	0,	2,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	8,	3,	2,	10,	8,	10,	9,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	11,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	11,	2,	8,	11,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	9,	0,	2,	3,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	11,	2,	1,	9,	11,	9,	8,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	10,	1,	11,	10,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	10,	1,	0,	8,	10,	8,	11,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	9,	0,	3,	11,	9,	11,	10,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	8,	10,	10,	8,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	7,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	3,	0,	7,	3,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	8,	4,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	1,	9,	4,	7,	1,	7,	3,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	8,	4,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	4,	7,	3,	0,	4,	1,	2,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	2,	10,	9,	0,	2,	8,	4,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	10,	9,	2,	9,	7,	2,	7,	3,	7,	9,	4,	-1,	-1,	-1,	-1	},
	{	8,	4,	7,	3,	11,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	4,	7,	11,	2,	4,	2,	0,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	0,	1,	8,	4,	7,	2,	3,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	7,	11,	9,	4,	11,	9,	11,	2,	9,	2,	1,	-1,	-1,	-1,	-1	},
	{	3,	10,	1,	3,	11,	10,	7,	8,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	11,	10,	1,	4,	11,	1,	0,	4,	7,	11,	4,	-1,	-1,	-1,	-1	},
	{	4,	7,	8,	9,	0,	11,	9,	11,	10,	11,	0,	3,	-1,	-1,	-1,	-1	},
	{	4,	7,	11,	4,	11,	9,	9,	11,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	5,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	5,	4,	0,	8,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	5,	4,	1,	5,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	5,	4,	8,	3,	5,	3,	1,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	9,	5,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	0,	8,	1,	2,	10,	4,	9,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	2,	10,	5,	4,	2,	4,	0,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	10,	5,	3,	2,	5,	3,	5,	4,	3,	4,	8,	-1,	-1,	-1,	-1	},
	{	9,	5,	4,	2,	3,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	11,	2,	0,	8,	11,	4,	9,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	5,	4,	0,	1,	5,	2,	3,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	1,	5,	2,	5,	8,	2,	8,	11,	4,	8,	5,	-1,	-1,	-1,	-1	},
	{	10,	3,	11,	10,	1,	3,	9,	5,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	9,	5,	0,	8,	1,	8,	10,	1,	8,	11,	10,	-1,	-1,	-1,	-1	},
	{	5,	4,	0,	5,	0,	11,	5,	11,	10,	11,	0,	3,	-1,	-1,	-1,	-1	},
	{	5,	4,	8,	5,	8,	10,	10,	8,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	7,	8,	5,	7,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	3,	0,	9,	5,	3,	5,	7,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	7,	8,	0,	1,	7,	1,	5,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	5,	3,	3,	5,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	7,	8,	9,	5,	7,	10,	1,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	1,	2,	9,	5,	0,	5,	3,	0,	5,	7,	3,	-1,	-1,	-1,	-1	},
	{	8,	0,	2,	8,	2,	5,	8,	5,	7,	10,	5,	2,	-1,	-1,	-1,	-1	},
	{	2,	10,	5,	2,	5,	3,	3,	5,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	9,	5,	7,	8,	9,	3,	11,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	5,	7,	9,	7,	2,	9,	2,	0,	2,	7,	11,	-1,	-1,	-1,	-1	},
	{	2,	3,	11,	0,	1,	8,	1,	7,	8,	1,	5,	7,	-1,	-1,	-1,	-1	},
	{	11,	2,	1,	11,	1,	7,	7,	1,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	5,	8,	8,	5,	7,	10,	1,	3,	10,	3,	11,	-1,	-1,	-1,	-1	},
	{	5,	7,	0,	5,	0,	9,	7,	11,	0,	1,	0,	10,	11,	10,	0,	-1	},
	{	11,	10,	0,	11,	0,	3,	10,	5,	0,	8,	0,	7,	5,	7,	0,	-1	},
	{	11,	10,	5,	7,	11,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	6,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	5,	10,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	0,	1,	5,	10,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	8,	3,	1,	9,	8,	5,	10,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	6,	5,	2,	6,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	6,	5,	1,	2,	6,	3,	0,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	6,	5,	9,	0,	6,	0,	2,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	9,	8,	5,	8,	2,	5,	2,	6,	3,	2,	8,	-1,	-1,	-1,	-1	},
	{	2,	3,	11,	10,	6,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	0,	8,	11,	2,	0,	10,	6,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	2,	3,	11,	5,	10,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	10,	6,	1,	9,	2,	9,	11,	2,	9,	8,	11,	-1,	-1,	-1,	-1	},
	{	6,	3,	11,	6,	5,	3,	5,	1,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	11,	0,	11,	5,	0,	5,	1,	5,	11,	6,	-1,	-1,	-1,	-1	},
	{	3,	11,	6,	0,	3,	6,	0,	6,	5,	0,	5,	9,	-1,	-1,	-1,	-1	},
	{	6,	5,	9,	6,	9,	11,	11,	9,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	10,	6,	4,	7,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	3,	0,	4,	7,	3,	6,	5,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	9,	0,	5,	10,	6,	8,	4,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	6,	5,	1,	9,	7,	1,	7,	3,	7,	9,	4,	-1,	-1,	-1,	-1	},
	{	6,	1,	2,	6,	5,	1,	4,	7,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	5,	5,	2,	6,	3,	0,	4,	3,	4,	7,	-1,	-1,	-1,	-1	},
	{	8,	4,	7,	9,	0,	5,	0,	6,	5,	0,	2,	6,	-1,	-1,	-1,	-1	},
	{	7,	3,	9,	7,	9,	4,	3,	2,	9,	5,	9,	6,	2,	6,	9,	-1	},
	{	3,	11,	2,	7,	8,	4,	10,	6,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	10,	6,	4,	7,	2,	4,	2,	0,	2,	7,	11,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	4,	7,	8,	2,	3,	11,	5,	10,	6,	-1,	-1,	-1,	-1	},
	{	9,	2,	1,	9,	11,	2,	9,	4,	11,	7,	11,	4,	5,	10,	6,	-1	},
	{	8,	4,	7,	3,	11,	5,	3,	5,	1,	5,	11,	6,	-1,	-1,	-1,	-1	},
	{	5,	1,	11,	5,	11,	6,	1,	0,	11,	7,	11,	4,	0,	4,	11,	-1	},
	{	0,	5,	9,	0,	6,	5,	0,	3,	6,	11,	6,	3,	8,	4,	7,	-1	},
	{	6,	5,	9,	6,	9,	11,	4,	7,	9,	7,	11,	9,	-1,	-1,	-1,	-1	},
	{	10,	4,	9,	6,	4,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	10,	6,	4,	9,	10,	0,	8,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	0,	1,	10,	6,	0,	6,	4,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	3,	1,	8,	1,	6,	8,	6,	4,	6,	1,	10,	-1,	-1,	-1,	-1	},
	{	1,	4,	9,	1,	2,	4,	2,	6,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	0,	8,	1,	2,	9,	2,	4,	9,	2,	6,	4,	-1,	-1,	-1,	-1	},
	{	0,	2,	4,	4,	2,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	3,	2,	8,	2,	4,	4,	2,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	4,	9,	10,	6,	4,	11,	2,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	2,	2,	8,	11,	4,	9,	10,	4,	10,	6,	-1,	-1,	-1,	-1	},
	{	3,	11,	2,	0,	1,	6,	0,	6,	4,	6,	1,	10,	-1,	-1,	-1,	-1	},
	{	6,	4,	1,	6,	1,	10,	4,	8,	1,	2,	1,	11,	8,	11,	1,	-1	},
	{	9,	6,	4,	9,	3,	6,	9,	1,	3,	11,	6,	3,	-1,	-1,	-1,	-1	},
	{	8,	11,	1,	8,	1,	0,	11,	6,	1,	9,	1,	4,	6,	4,	1,	-1	},
	{	3,	11,	6,	3,	6,	0,	0,	6,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	6,	4,	8,	11,	6,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	10,	6,	7,	8,	10,	8,	9,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	7,	3,	0,	10,	7,	0,	9,	10,	6,	7,	10,	-1,	-1,	-1,	-1	},
	{	10,	6,	7,	1,	10,	7,	1,	7,	8,	1,	8,	0,	-1,	-1,	-1,	-1	},
	{	10,	6,	7,	10,	7,	1,	1,	7,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	6,	1,	6,	8,	1,	8,	9,	8,	6,	7,	-1,	-1,	-1,	-1	},
	{	2,	6,	9,	2,	9,	1,	6,	7,	9,	0,	9,	3,	7,	3,	9,	-1	},
	{	7,	8,	0,	7,	0,	6,	6,	0,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	3,	2,	6,	7,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	3,	11,	10,	6,	8,	10,	8,	9,	8,	6,	7,	-1,	-1,	-1,	-1	},
	{	2,	0,	7,	2,	7,	11,	0,	9,	7,	6,	7,	10,	9,	10,	7,	-1	},
	{	1,	8,	0,	1,	7,	8,	1,	10,	7,	6,	7,	10,	2,	3,	11,	-1	},
	{	11,	2,	1,	11,	1,	7,	10,	6,	1,	6,	7,	1,	-1,	-1,	-1,	-1	},
	{	8,	9,	6,	8,	6,	7,	9,	1,	6,	11,	6,	3,	1,	3,	6,	-1	},
	{	0,	9,	1,	11,	6,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	8,	0,	7,	0,	6,	3,	11,	0,	11,	6,	0,	-1,	-1,	-1,	-1	},
	{	7,	11,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	6,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	0,	8,	11,	7,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	11,	7,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	1,	9,	8,	3,	1,	11,	7,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	1,	2,	6,	11,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	3,	0,	8,	6,	11,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	9,	0,	2,	10,	9,	6,	11,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	6,	11,	7,	2,	10,	3,	10,	8,	3,	10,	9,	8,	-1,	-1,	-1,	-1	},
	{	7,	2,	3,	6,	2,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	7,	0,	8,	7,	6,	0,	6,	2,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	7,	6,	2,	3,	7,	0,	1,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	6,	2,	1,	8,	6,	1,	9,	8,	8,	7,	6,	-1,	-1,	-1,	-1	},
	{	10,	7,	6,	10,	1,	7,	1,	3,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	7,	6,	1,	7,	10,	1,	8,	7,	1,	0,	8,	-1,	-1,	-1,	-1	},
	{	0,	3,	7,	0,	7,	10,	0,	10,	9,	6,	10,	7,	-1,	-1,	-1,	-1	},
	{	7,	6,	10,	7,	10,	8,	8,	10,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	6,	8,	4,	11,	8,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	6,	11,	3,	0,	6,	0,	4,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	6,	11,	8,	4,	6,	9,	0,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	4,	6,	9,	6,	3,	9,	3,	1,	11,	3,	6,	-1,	-1,	-1,	-1	},
	{	6,	8,	4,	6,	11,	8,	2,	10,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	3,	0,	11,	0,	6,	11,	0,	4,	6,	-1,	-1,	-1,	-1	},
	{	4,	11,	8,	4,	6,	11,	0,	2,	9,	2,	10,	9,	-1,	-1,	-1,	-1	},
	{	10,	9,	3,	10,	3,	2,	9,	4,	3,	11,	3,	6,	4,	6,	3,	-1	},
	{	8,	2,	3,	8,	4,	2,	4,	6,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	4,	2,	4,	6,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	9,	0,	2,	3,	4,	2,	4,	6,	4,	3,	8,	-1,	-1,	-1,	-1	},
	{	1,	9,	4,	1,	4,	2,	2,	4,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	1,	3,	8,	6,	1,	8,	4,	6,	6,	10,	1,	-1,	-1,	-1,	-1	},
	{	10,	1,	0,	10,	0,	6,	6,	0,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	6,	3,	4,	3,	8,	6,	10,	3,	0,	3,	9,	10,	9,	3,	-1	},
	{	10,	9,	4,	6,	10,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	9,	5,	7,	6,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	4,	9,	5,	11,	7,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	0,	1,	5,	4,	0,	7,	6,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	7,	6,	8,	3,	4,	3,	5,	4,	3,	1,	5,	-1,	-1,	-1,	-1	},
	{	9,	5,	4,	10,	1,	2,	7,	6,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	6,	11,	7,	1,	2,	10,	0,	8,	3,	4,	9,	5,	-1,	-1,	-1,	-1	},
	{	7,	6,	11,	5,	4,	10,	4,	2,	10,	4,	0,	2,	-1,	-1,	-1,	-1	},
	{	3,	4,	8,	3,	5,	4,	3,	2,	5,	10,	5,	2,	11,	7,	6,	-1	},
	{	7,	2,	3,	7,	6,	2,	5,	4,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	5,	4,	0,	8,	6,	0,	6,	2,	6,	8,	7,	-1,	-1,	-1,	-1	},
	{	3,	6,	2,	3,	7,	6,	1,	5,	0,	5,	4,	0,	-1,	-1,	-1,	-1	},
	{	6,	2,	8,	6,	8,	7,	2,	1,	8,	4,	8,	5,	1,	5,	8,	-1	},
	{	9,	5,	4,	10,	1,	6,	1,	7,	6,	1,	3,	7,	-1,	-1,	-1,	-1	},
	{	1,	6,	10,	1,	7,	6,	1,	0,	7,	8,	7,	0,	9,	5,	4,	-1	},
	{	4,	0,	10,	4,	10,	5,	0,	3,	10,	6,	10,	7,	3,	7,	10,	-1	},
	{	7,	6,	10,	7,	10,	8,	5,	4,	10,	4,	8,	10,	-1,	-1,	-1,	-1	},
	{	6,	9,	5,	6,	11,	9,	11,	8,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	6,	11,	0,	6,	3,	0,	5,	6,	0,	9,	5,	-1,	-1,	-1,	-1	},
	{	0,	11,	8,	0,	5,	11,	0,	1,	5,	5,	6,	11,	-1,	-1,	-1,	-1	},
	{	6,	11,	3,	6,	3,	5,	5,	3,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	10,	9,	5,	11,	9,	11,	8,	11,	5,	6,	-1,	-1,	-1,	-1	},
	{	0,	11,	3,	0,	6,	11,	0,	9,	6,	5,	6,	9,	1,	2,	10,	-1	},
	{	11,	8,	5,	11,	5,	6,	8,	0,	5,	10,	5,	2,	0,	2,	5,	-1	},
	{	6,	11,	3,	6,	3,	5,	2,	10,	3,	10,	5,	3,	-1,	-1,	-1,	-1	},
	{	5,	8,	9,	5,	2,	8,	5,	6,	2,	3,	8,	2,	-1,	-1,	-1,	-1	},
	{	9,	5,	6,	9,	6,	0,	0,	6,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	5,	8,	1,	8,	0,	5,	6,	8,	3,	8,	2,	6,	2,	8,	-1	},
	{	1,	5,	6,	2,	1,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	3,	6,	1,	6,	10,	3,	8,	6,	5,	6,	9,	8,	9,	6,	-1	},
	{	10,	1,	0,	10,	0,	6,	9,	5,	0,	5,	6,	0,	-1,	-1,	-1,	-1	},
	{	0,	3,	8,	5,	6,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	5,	6,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	5,	10,	7,	5,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	5,	10,	11,	7,	5,	8,	3,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	11,	7,	5,	10,	11,	1,	9,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	10,	7,	5,	10,	11,	7,	9,	8,	1,	8,	3,	1,	-1,	-1,	-1,	-1	},
	{	11,	1,	2,	11,	7,	1,	7,	5,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	1,	2,	7,	1,	7,	5,	7,	2,	11,	-1,	-1,	-1,	-1	},
	{	9,	7,	5,	9,	2,	7,	9,	0,	2,	2,	11,	7,	-1,	-1,	-1,	-1	},
	{	7,	5,	2,	7,	2,	11,	5,	9,	2,	3,	2,	8,	9,	8,	2,	-1	},
	{	2,	5,	10,	2,	3,	5,	3,	7,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	2,	0,	8,	5,	2,	8,	7,	5,	10,	2,	5,	-1,	-1,	-1,	-1	},
	{	9,	0,	1,	5,	10,	3,	5,	3,	7,	3,	10,	2,	-1,	-1,	-1,	-1	},
	{	9,	8,	2,	9,	2,	1,	8,	7,	2,	10,	2,	5,	7,	5,	2,	-1	},
	{	1,	3,	5,	3,	7,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	7,	0,	7,	1,	1,	7,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	0,	3,	9,	3,	5,	5,	3,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	8,	7,	5,	9,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	8,	4,	5,	10,	8,	10,	11,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	5,	0,	4,	5,	11,	0,	5,	10,	11,	11,	3,	0,	-1,	-1,	-1,	-1	},
	{	0,	1,	9,	8,	4,	10,	8,	10,	11,	10,	4,	5,	-1,	-1,	-1,	-1	},
	{	10,	11,	4,	10,	4,	5,	11,	3,	4,	9,	4,	1,	3,	1,	4,	-1	},
	{	2,	5,	1,	2,	8,	5,	2,	11,	8,	4,	5,	8,	-1,	-1,	-1,	-1	},
	{	0,	4,	11,	0,	11,	3,	4,	5,	11,	2,	11,	1,	5,	1,	11,	-1	},
	{	0,	2,	5,	0,	5,	9,	2,	11,	5,	4,	5,	8,	11,	8,	5,	-1	},
	{	9,	4,	5,	2,	11,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	5,	10,	3,	5,	2,	3,	4,	5,	3,	8,	4,	-1,	-1,	-1,	-1	},
	{	5,	10,	2,	5,	2,	4,	4,	2,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	10,	2,	3,	5,	10,	3,	8,	5,	4,	5,	8,	0,	1,	9,	-1	},
	{	5,	10,	2,	5,	2,	4,	1,	9,	2,	9,	4,	2,	-1,	-1,	-1,	-1	},
	{	8,	4,	5,	8,	5,	3,	3,	5,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	4,	5,	1,	0,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	8,	4,	5,	8,	5,	3,	9,	0,	5,	0,	3,	5,	-1,	-1,	-1,	-1	},
	{	9,	4,	5,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	11,	7,	4,	9,	11,	9,	10,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	8,	3,	4,	9,	7,	9,	11,	7,	9,	10,	11,	-1,	-1,	-1,	-1	},
	{	1,	10,	11,	1,	11,	4,	1,	4,	0,	7,	4,	11,	-1,	-1,	-1,	-1	},
	{	3,	1,	4,	3,	4,	8,	1,	10,	4,	7,	4,	11,	10,	11,	4,	-1	},
	{	4,	11,	7,	9,	11,	4,	9,	2,	11,	9,	1,	2,	-1,	-1,	-1,	-1	},
	{	9,	7,	4,	9,	11,	7,	9,	1,	11,	2,	11,	1,	0,	8,	3,	-1	},
	{	11,	7,	4,	11,	4,	2,	2,	4,	0,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	11,	7,	4,	11,	4,	2,	8,	3,	4,	3,	2,	4,	-1,	-1,	-1,	-1	},
	{	2,	9,	10,	2,	7,	9,	2,	3,	7,	7,	4,	9,	-1,	-1,	-1,	-1	},
	{	9,	10,	7,	9,	7,	4,	10,	2,	7,	8,	7,	0,	2,	0,	7,	-1	},
	{	3,	7,	10,	3,	10,	2,	7,	4,	10,	1,	10,	0,	4,	0,	10,	-1	},
	{	1,	10,	2,	8,	7,	4,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	9,	1,	4,	1,	7,	7,	1,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	9,	1,	4,	1,	7,	0,	8,	1,	8,	7,	1,	-1,	-1,	-1,	-1	},
	{	4,	0,	3,	7,	4,	3,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	4,	8,	7,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	10,	8,	10,	11,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	0,	9,	3,	9,	11,	11,	9,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	1,	10,	0,	10,	8,	8,	10,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	1,	10,	11,	3,	10,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	2,	11,	1,	11,	9,	9,	11,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	0,	9,	3,	9,	11,	1,	2,	9,	2,	11,	9,	-1,	-1,	-1,	-1	},
	{	0,	2,	11,	8,	0,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	3,	2,	11,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	3,	8,	2,	8,	10,	10,	8,	9,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	9,	10,	2,	0,	9,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	2,	3,	8,	2,	8,	10,	0,	1,	8,	1,	10,	8,	-1,	-1,	-1,	-1	},
	{	1,	10,	2,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	1,	3,	8,	9,	1,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	9,	1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	0,	3,	8,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	},
	{	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1,	-1	}
};


static inline point3D binInterp(int x1, int y1, int z1, int v1, int x2, int y2, int z2, int v2) {
//	if (v1 == v2)
//		return point3D((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0);
//	else
//	if(v1)
//		return point3D(x1, y1, z1);
//	else
//		return point3D(x2, y2, z2);
	return point3D((x1+x2)/2.0, (y1+y2)/2.0, (z1+z2)/2.0);
}

static inline void polygonise(uint16_t i, uint16_t j, uint16_t k,
		bool const * inside_surface, uint16_t w, uint16_t h,
		vector<faceinfo> & faces, vector<point3D> & vertices,
		unordered_map<point3D, uint32_t> & vertex_map) {
	point3D vertlist[12];
	/*
		Determine the index into the edge table which
		tells us which vertices are inside of the surface
	*/
	int cubeindex = 0;
	if (inside_surface[( i   *w+j)  *h+k  ]) cubeindex |= 1;
	if (inside_surface[((i+1)*w+j)  *h+k  ]) cubeindex |= 2;
	if (inside_surface[((i+1)*w+j+1)*h+k  ]) cubeindex |= 4;
	if (inside_surface[( i   *w+j+1)*h+k  ]) cubeindex |= 8;
	if (inside_surface[( i   *w+j)  *h+k+1]) cubeindex |= 16;
	if (inside_surface[((i+1)*w+j)  *h+k+1]) cubeindex |= 32;
	if (inside_surface[((i+1)*w+j+1)*h+k+1]) cubeindex |= 64;
	if (inside_surface[( i   *w+j+1)*h+k+1]) cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		return;

	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
		vertlist[0] = point3D(i + 0.5, j, k);
	if (edgeTable[cubeindex] & 2)
		vertlist[1] = point3D(i + 1, j + 0.5, k);
	if (edgeTable[cubeindex] & 4)
		vertlist[2] = point3D(i + 0.5, j + 1, k);
	if (edgeTable[cubeindex] & 8)
		vertlist[3] = point3D(i, j + 0.5, k);
	if (edgeTable[cubeindex] & 16)
		vertlist[4] = point3D(i + 0.5, j, k + 1);
	if (edgeTable[cubeindex] & 32)
		vertlist[5] = point3D(i + 1, j + 0.5, k + 1);
	if (edgeTable[cubeindex] & 64)
		vertlist[6] = point3D(i + 0.5, j + 1, k + 1);
	if (edgeTable[cubeindex] & 128)
		vertlist[7] = point3D(i, j + 0.5, k + 1);
	if (edgeTable[cubeindex] & 256)
		vertlist[8] = point3D(i, j, k + 0.5);
	if (edgeTable[cubeindex] & 512)
		vertlist[9] = point3D(i + 1, j, k + 0.5);
	if (edgeTable[cubeindex] & 1024)
		vertlist[10] = point3D(i + 1, j + 1, k + 0.5);
	if (edgeTable[cubeindex] & 2048)
		vertlist[11] = point3D(i, j + 1, k + 0.5);

	point3D p1, p2, p3;
	uint32_t id1, id2, id3;
	for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {
		p1 = vertlist[triTable[cubeindex][i  ]];
		p2 = vertlist[triTable[cubeindex][i+1]];
		p3 = vertlist[triTable[cubeindex][i+2]];

		if (vertex_map.find(p1) == vertex_map.end()) {
			vertices.push_back(p1);
			id1 = vertices.size() - 1;
			vertex_map[p1] = id1;
		}
		else {
			id1 = vertex_map[p1];
		}
		if (vertex_map.find(p2) == vertex_map.end()) {
			vertices.push_back(p2);
			id2 = vertices.size() - 1;
			vertex_map[p2] = id2;
		}
		else {
			id2 = vertex_map[p2];
		}
		if (vertex_map.find(p3) == vertex_map.end()) {
			vertices.push_back(p3);
			id3 = vertices.size() - 1;
			vertex_map[p3] = id3;
		}
		else {
			id3 = vertex_map[p3];
		}
		faces.push_back({id1, id2, id3});
	}
}

static inline void polygonise(uint16_t i, uint16_t j, uint16_t k,
		bool const * inside_surface, uint16_t w, uint16_t h,
		vector<triangle> & triangles) {
	point3D vertlist[12];
	/*
		Determine the index into the edge table which
		tells us which vertices are inside of the surface
	*/
	int cubeindex = 0;
	if (inside_surface[( i   *w+j)  *h+k  ]) cubeindex |= 1;
	if (inside_surface[((i+1)*w+j)  *h+k  ]) cubeindex |= 2;
	if (inside_surface[((i+1)*w+j+1)*h+k  ]) cubeindex |= 4;
	if (inside_surface[( i   *w+j+1)*h+k  ]) cubeindex |= 8;
	if (inside_surface[( i   *w+j)  *h+k+1]) cubeindex |= 16;
	if (inside_surface[((i+1)*w+j)  *h+k+1]) cubeindex |= 32;
	if (inside_surface[((i+1)*w+j+1)*h+k+1]) cubeindex |= 64;
	if (inside_surface[( i   *w+j+1)*h+k+1]) cubeindex |= 128;

	/* Cube is entirely in/out of the surface */
	if (edgeTable[cubeindex] == 0)
		return;

	/* Find the vertices where the surface intersects the cube */
	if (edgeTable[cubeindex] & 1)
		vertlist[0] = point3D(i + 0.5, j, k);
	if (edgeTable[cubeindex] & 2)
		vertlist[1] = point3D(i + 1, j + 0.5, k);
	if (edgeTable[cubeindex] & 4)
		vertlist[2] = point3D(i + 0.5, j + 1, k);
	if (edgeTable[cubeindex] & 8)
		vertlist[3] = point3D(i, j + 0.5, k);
	if (edgeTable[cubeindex] & 16)
		vertlist[4] = point3D(i + 0.5, j, k + 1);
	if (edgeTable[cubeindex] & 32)
		vertlist[5] = point3D(i + 1, j + 0.5, k + 1);
	if (edgeTable[cubeindex] & 64)
		vertlist[6] = point3D(i + 0.5, j + 1, k + 1);
	if (edgeTable[cubeindex] & 128)
		vertlist[7] = point3D(i, j + 0.5, k + 1);
	if (edgeTable[cubeindex] & 256)
		vertlist[8] = point3D(i, j, k + 0.5);
	if (edgeTable[cubeindex] & 512)
		vertlist[9] = point3D(i + 1, j, k + 0.5);
	if (edgeTable[cubeindex] & 1024)
		vertlist[10] = point3D(i + 1, j + 1, k + 0.5);
	if (edgeTable[cubeindex] & 2048)
		vertlist[11] = point3D(i, j + 1, k + 0.5);

	point3D p1, p2, p3;
	uint32_t id1, id2, id3;
	for (int i = 0; triTable[cubeindex][i] != -1; i += 3) {
		p1 = vertlist[triTable[cubeindex][i  ]];
		p2 = vertlist[triTable[cubeindex][i+1]];
		p3 = vertlist[triTable[cubeindex][i+2]];
		triangles.push_back({p1, p2, p3});
	}
}

static void marchingCubes(bool const * cpk_model, uint16_t pl, uint16_t pw, uint16_t ph,
		vector<faceinfo> & faceList, vector<point3D> & vertList) {
	// Polygonise the grid
	unordered_map<point3D, uint32_t> vertex_map;

	vector<triangle> triangle_vector;
	for (uint16_t i = 0; i < pl - 1; ++i) {
		for (uint16_t j = 0; j < pw - 1; ++j) {
			for (uint16_t k = 0; k < ph - 1; ++k) {
				if (cpk_model[(i  * pw + j  ) * ph + k]
				 || cpk_model[(i  * pw + j+1) * ph + k]
				 || cpk_model[(i  * pw + j  ) * ph + k+1]
				 || cpk_model[(i  * pw + j+1) * ph + k+1]
				 || cpk_model[((i+1)* pw + j  ) * ph + k]
				 || cpk_model[((i+1)* pw + j+1) * ph + k]
				 || cpk_model[((i+1)* pw + j  ) * ph + k+1]
				 || cpk_model[((i+1)* pw + j+1) * ph + k+1])
					polygonise(i, j, k, cpk_model, pw, ph, faceList, vertList, vertex_map);
			}
		}
	}
}



static void marchingCubesParallel(bool const * cpk_model, uint16_t pl, uint16_t pw, uint16_t ph,
		vector<faceinfo> & faceList, vector<point3D> & vertList) {
	const int nthreads = omp_get_max_threads();

	vector<vector<triangle>> triangles(nthreads);
#pragma omp parallel
{
	const int ithread = omp_get_thread_num();
	#pragma omp for schedule(runtime)
	for (uint16_t i = 0; i < pl - 1; ++i) {
		for (uint16_t j = 0; j < pw - 1; ++j) {
			for (uint16_t k = 0; k < ph - 1; ++k) {
				if (cpk_model[(i  * pw + j  ) * ph + k]
				 || cpk_model[(i  * pw + j+1) * ph + k]
				 || cpk_model[(i  * pw + j  ) * ph + k+1]
				 || cpk_model[(i  * pw + j+1) * ph + k+1]
				 || cpk_model[((i+1)* pw + j  ) * ph + k]
				 || cpk_model[((i+1)* pw + j+1) * ph + k]
				 || cpk_model[((i+1)* pw + j  ) * ph + k+1]
				 || cpk_model[((i+1)* pw + j+1) * ph + k+1])
					polygonise(i, j, k, cpk_model, pw, ph, triangles[ithread]);
			}
		}
	}
}

	unordered_map<point3D, uint32_t> vertMap;
	for (size_t ii = 0; ii < nthreads; ++ii) {
		for (auto const & t : triangles[ii]) {
			uint32_t id1, id2, id3;
			if (vertMap.find(t.p1) == vertMap.end()) {
				vertList.push_back(t.p1);
				id1 = vertList.size() - 1;
				vertMap[t.p1] = id1;
			} else {
				id1 = vertMap[t.p1];
			}
			if (vertMap.find(t.p2) == vertMap.end()) {
				vertList.push_back(t.p2);
				id2 = vertList.size() - 1;
				vertMap[t.p2] = id2;
			} else {
				id2 = vertMap[t.p2];
			}
			if (vertMap.find(t.p3) == vertMap.end()) {
				vertList.push_back(t.p3);
				id3 = vertList.size() - 1;
				vertMap[t.p3] = id3;
			} else {
				id3 = vertMap[t.p3];
			}
			faceList.push_back({id1, id2, id3});
		}
	}
}



static void marchingCubesParallel2(bool const * cpk_model, vector<atom> const * atoms,
		uint16_t pl, uint16_t pw, uint16_t ph, point3D const & ptran, float probeRadius, float resolution,
		vector<faceinfo> & faceList, vector<point3D> & vertList) {
	const int nthreads = omp_get_max_threads();

	vector<vector<faceinfo>> faces(nthreads);
	vector<vector<point3D>> vertices(nthreads);
	vector<unordered_map<point3D, uint32_t>> vertex_map(nthreads);
#pragma omp parallel
{
	const int ithread = omp_get_thread_num();
	#pragma omp for schedule(runtime)
	for (size_t i = 0; i < atoms->size(); ++i) {
		/* Translate and discretize the coordinates */
		int cx = static_cast<int>(round(((*atoms)[i].x + ptran.x) * resolution));
		int cy = static_cast<int>(round(((*atoms)[i].y + ptran.y) * resolution));
		int cz = static_cast<int>(round(((*atoms)[i].z + ptran.z) * resolution));
		int c  = static_cast<int>(round(((*atoms)[i].radius + probeRadius) * resolution + 0.5));

		for (int ii = cx - c; ii <= cx + c; ++ii) {
			for (int jj = cy - c; jj <= cy + c; ++jj) {
				for (int kk = cz - c; kk <= cz + c; ++kk) {
					polygonise(ii, jj, kk, cpk_model, pw, ph, faces[ithread], vertices[ithread], vertex_map[ithread]);
				}
			}
		}
	}
}
	vertList = vertices[0];
	faceList = faces[0];
	unordered_map<point3D, uint32_t> vertMap = vertex_map[0];
	for (size_t ii = 1; ii < nthreads; ++ii) {
		for (auto const & f : faces[ii]) {
			uint32_t id1, id2, id3;
			if (vertMap.count(vertices[ii][f.a]) == 0) {
				vertList.push_back(vertices[ii][f.a]);
				id1 = vertList.size() - 1;
				vertMap[vertices[ii][f.a]] = id1;
			} else {
				id1 = vertMap[vertices[ii][f.a]];
			}
			if (vertMap.count(vertices[ii][f.b]) == 0) {
				vertList.push_back(vertices[ii][f.b]);
				id2 = vertList.size() - 1;
				vertMap[vertices[ii][f.b]] = id2;
			} else {
				id2 = vertMap[vertices[ii][f.b]];
			}
			if (vertMap.count(vertices[ii][f.c]) == 0) {
				vertList.push_back(vertices[ii][f.c]);
				id3 = vertList.size() - 1;
				vertMap[vertices[ii][f.c]] = id3;
			} else {
				id3 = vertMap[vertices[ii][f.c]];
			}
			faceList.push_back({id1, id2, id3});
		}
	}
}

#endif /* MARCHINGCUBES_H_ */