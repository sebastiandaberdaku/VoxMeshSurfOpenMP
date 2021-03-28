/*
 * DrawSphere.h
 *
 *  Created on: 02/feb/2015
 *      Author: sebastian
 */

#ifndef SPHERE_DRAWSPHERE_H_
#define SPHERE_DRAWSPHERE_H_


/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the XY plane). This method gets an initial
 * error value as an additional parameter, in order to compensate for the rounding error of
 * the center coordinates and radius.
 *
 * If no initial error is propagated the sphere would raster to a cylinder near its equator,
 * as the same circle is drawn several times, where in reality the arc should vary slightly
 * even if the x=0 or y=0 points would stay the same.
 *
 * Using the accumulated error from each point as an input to the next circle solves the issue.
 */
inline static void DrawCircleXY(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y0, int z, int radius, int error0) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		grid[((x0 + x) * pwidth + (y0 + y)) * pheight + z]= true;
		grid[((x0 + x) * pwidth + (y0 - y)) * pheight + z]= true;
		grid[((x0 - x) * pwidth + (y0 + y)) * pheight + z]= true;
		grid[((x0 - x) * pwidth + (y0 - y)) * pheight + z]= true;

		grid[((x0 + y) * pwidth + (y0 + x)) * pheight + z]= true;
		grid[((x0 + y) * pwidth + (y0 - x)) * pheight + z]= true;
		grid[((x0 - y) * pwidth + (y0 + x)) * pheight + z]= true;
		grid[((x0 - y) * pwidth + (y0 - x)) * pheight + z]= true;

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};

inline static void DrawDiskXY(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y0, int z, int radius, int error0) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		for (int i = y; i >= x; --i) {
			grid[((x0 + x) * pwidth + (y0 + i)) * pheight + z]= true;
			grid[((x0 - x) * pwidth + (y0 + i)) * pheight + z]= true;
			grid[((x0 - i) * pwidth + (y0 + x)) * pheight + z]= true;
			grid[((x0 + i) * pwidth + (y0 + x)) * pheight + z]= true;

			grid[((x0 + x) * pwidth + (y0 - i)) * pheight + z]= true;
			grid[((x0 - x) * pwidth + (y0 - i)) * pheight + z]= true;
			grid[((x0 - i) * pwidth + (y0 - x)) * pheight + z]= true;
			grid[((x0 + i) * pwidth + (y0 - x)) * pheight + z]= true;
		}

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};

inline static void DrawDiskXY(int z, int radius, int error0, vector<voxel_offset> & ball_offsets) {
	int y = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (y >= x) {
		for (int i = y; i >= x; --i) {
			ball_offsets.push_back({ x, i, z});
			ball_offsets.push_back({-x, i, z});
			ball_offsets.push_back({ i, x, z});
			ball_offsets.push_back({-i, x, z});

			ball_offsets.push_back({ x,-i, z});
			ball_offsets.push_back({-x,-i, z});
			ball_offsets.push_back({ i,-x, z});
			ball_offsets.push_back({-i,-x, z});
		}

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--y;
			radiusError += 2 * (x - y) + 1;
		}
	}
};

/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the YZ plane).
 */
inline static void DrawCircleYZ(bool * grid, uint16_t pwidth, uint16_t pheight, int x, int y0, int z0, int radius, int error0) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		grid[(x * pwidth + (y0 + y)) * pheight + z0 + z]= true;
		grid[(x * pwidth + (y0 + y)) * pheight + z0 - z]= true;
		grid[(x * pwidth + (y0 - y)) * pheight + z0 + z]= true;
		grid[(x * pwidth + (y0 - y)) * pheight + z0 - z]= true;

		grid[(x * pwidth + (y0 + z)) * pheight + z0 + y]= true;
		grid[(x * pwidth + (y0 + z)) * pheight + z0 - y]= true;
		grid[(x * pwidth + (y0 - z)) * pheight + z0 + y]= true;
		grid[(x * pwidth + (y0 - z)) * pheight + z0 - y]= true;

		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};

inline static void DrawCircleYZf(bool * grid, uint16_t pwidth, uint16_t pheight, float x, float y0, float z0, float radius, float error0) {
	float z = radius, y = 0;
	float radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		grid[static_cast<int>((x * pwidth + (y0 + y)) * pheight + z0 + z)]= true;
		grid[static_cast<int>((x * pwidth + (y0 + y)) * pheight + z0 - z)]= true;
		grid[static_cast<int>((x * pwidth + (y0 - y)) * pheight + z0 + z)]= true;
		grid[static_cast<int>((x * pwidth + (y0 - y)) * pheight + z0 - z)]= true;

		grid[static_cast<int>((x * pwidth + (y0 + z)) * pheight + z0 + y)]= true;
		grid[static_cast<int>((x * pwidth + (y0 + z)) * pheight + z0 - y)]= true;
		grid[static_cast<int>((x * pwidth + (y0 - z)) * pheight + z0 + y)]= true;
		grid[static_cast<int>((x * pwidth + (y0 - z)) * pheight + z0 - y)]= true;

		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};

inline static void DrawCircleYZ(int x, int radius, int error0, vector<voxel_offset> & ball_offsets) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		ball_offsets.push_back({ x, y, z});
		ball_offsets.push_back({ x, y,-z});
		ball_offsets.push_back({ x,-y, z});
		ball_offsets.push_back({ x,-y,-z});

		ball_offsets.push_back({ x, z, y});
		ball_offsets.push_back({ x, z,-y});
		ball_offsets.push_back({ x,-z, y});
		ball_offsets.push_back({ x,-z,-y});

		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};

/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the YZ plane).
 */
inline static void DrawDiskYZ(bool * grid, uint16_t pwidth, uint16_t pheight, int x, int y0, int z0, int radius, int error0) {
	int z = radius, y = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= y) {
		for (int i = z; i >= y; --i) {
			grid[(x * pwidth + (y0 + i)) * pheight + z0 + y]= true;
			grid[(x * pwidth + (y0 + i)) * pheight + z0 - y]= true;
			grid[(x * pwidth + (y0 - i)) * pheight + z0 + y]= true;
			grid[(x * pwidth + (y0 - i)) * pheight + z0 - y]= true;

			grid[(x * pwidth + (y0 + y)) * pheight + z0 + i]= true;
			grid[(x * pwidth + (y0 + y)) * pheight + z0 - i]= true;
			grid[(x * pwidth + (y0 - y)) * pheight + z0 + i]= true;
			grid[(x * pwidth + (y0 - y)) * pheight + z0 - i]= true;
		}
		++y;
		if (radiusError < 0) {
			radiusError += 2 * y + 1;
		} else {
			--z;
			radiusError += 2 * (y - z) + 1;
		}
	}
};



/**
 * Implementation of the Midpoint Circle Algorithm which determines the points needed for
 * drawing a circle on a given plane (parallel to the XZ plane).
 */
inline static void DrawCircleXZ(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y, int z0, int radius, int error0) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		grid[((x0 + x) * pwidth + y) * pheight + z0 + z]= true;
		grid[((x0 + x) * pwidth + y) * pheight + z0 - z]= true;
		grid[((x0 - x) * pwidth + y) * pheight + z0 + z]= true;
		grid[((x0 - x) * pwidth + y) * pheight + z0 - z]= true;

		grid[((x0 + z) * pwidth + y) * pheight + z0 + x]= true;
		grid[((x0 + z) * pwidth + y) * pheight + z0 - x]= true;
		grid[((x0 - z) * pwidth + y) * pheight + z0 + x]= true;
		grid[((x0 - z) * pwidth + y) * pheight + z0 - x]= true;

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};

inline static void DrawCircleXZf(bool * grid, uint16_t pwidth, uint16_t pheight, float x0, float y, float z0, float radius, float error0) {
	float z = radius, x = 0;
	float radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		grid[static_cast<int>(((x0 + x) * pwidth + y) * pheight + z0 + z)]= true;
		grid[static_cast<int>(((x0 + x) * pwidth + y) * pheight + z0 - z)]= true;
		grid[static_cast<int>(((x0 - x) * pwidth + y) * pheight + z0 + z)]= true;
		grid[static_cast<int>(((x0 - x) * pwidth + y) * pheight + z0 - z)]= true;

		grid[static_cast<int>(((x0 + z) * pwidth + y) * pheight + z0 + x)]= true;
		grid[static_cast<int>(((x0 + z) * pwidth + y) * pheight + z0 - x)]= true;
		grid[static_cast<int>(((x0 - z) * pwidth + y) * pheight + z0 + x)]= true;
		grid[static_cast<int>(((x0 - z) * pwidth + y) * pheight + z0 - x)]= true;

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};

inline static void DrawCircleXZ(int y, int radius, int error0, vector<voxel_offset> & ball_offsets) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		ball_offsets.push_back({ x, y, z});
		ball_offsets.push_back({ x, y,-z});
		ball_offsets.push_back({-x, y, z});
		ball_offsets.push_back({-x, y,-z});

		ball_offsets.push_back({ z, y, x});
		ball_offsets.push_back({ z, y,-x});
		ball_offsets.push_back({-z, y, x});
		ball_offsets.push_back({-z, y,-x});

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};

inline static void DrawDiskXZ(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y, int z0, int radius, int error0) {
	int z = radius, x = 0;
	int radiusError = error0; // Initial error state passed in, NOT 1-x

	while (z >= x) {
		for (int i = z0 - z; i <= z0 + z; ++i) {
			grid[((x0 + x) * pwidth + y) * pheight + i]= true;
			grid[((x0 - x) * pwidth + y) * pheight + i]= true;
		}
		for (int i = z0 - x; i <= z0 + x; ++i) {
			grid[((x0 + z) * pwidth + y) * pheight + i]= true;
			grid[((x0 - z) * pwidth + y) * pheight + i]= true;
		}

		++x;
		if (radiusError < 0) {
			radiusError += 2 * x + 1;
		} else {
			--z;
			radiusError += 2 * (x - z) + 1;
		}
	}
};
/**
 * Adapted Midpoint Circle Algorithm for drawing spheres in 3D.
 * @param grid		target voxelGrid
 * @param x0		x coordinate of the center of the sphere
 * @param y0		y coordinate of the center of the sphere
 * @param z0		z coordinate of the center of the sphere
 * @param radius	radius of the sphere
 */

inline static void DrawSphere(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y0, int z0, int radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawCircleXY(grid, pwidth, pheight, x0, y0, z0 - t, r, radiusError);
		DrawCircleXY(grid, pwidth, pheight, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 - t, z0, r, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};

inline static void DrawBall(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y0, int z0, int radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 - r, t, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 + r, t, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 - t, r, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 - t, z0, r, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};

inline static void DrawBalln(bool * grid, uint16_t pwidth, uint16_t pheight, int x0, int y0, int z0, int radius) {
	for (int i = -radius; i <= radius; ++i) {
		for (int j = -radius; j <= radius; ++j) {
			for (int k = -radius; k <= radius; ++k) {
				if (i*i + j*j + k*k <= radius * radius) {
					grid[((x0 + i) * pwidth + (y0 + j)) * pheight + (z0 + k)]= true;
				}
			}
		}
	}
};


//inline static void DrawDiskXYf(bool * grid, uint16_t pwidth, uint16_t pheight, float x0, float y0, float z, float radius, float error0) {
//	float y = radius, x = 0;
//	float radiusError = error0; // Initial error state passed in, NOT 1-x
//
//	while (y >= x) {
//		for (float i = y; i >= x; --i) {
//			grid[static_cast<int>(((x0 + x) * pwidth + (y0 + i)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 - x) * pwidth + (y0 + i)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 - i) * pwidth + (y0 + x)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 + i) * pwidth + (y0 + x)) * pheight + z)]= true;
//
//			grid[static_cast<int>(((x0 + x) * pwidth + (y0 - i)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 - x) * pwidth + (y0 - i)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 - i) * pwidth + (y0 - x)) * pheight + z)]= true;
//			grid[static_cast<int>(((x0 + i) * pwidth + (y0 - x)) * pheight + z)]= true;
//		}
//
//		++x;
//		if (radiusError < 0) {
//			radiusError += 2 * x + 1;
//		} else {
//			--y;
//			radiusError += 2 * (x - y) + 1;
//		}
//	}
//};

inline static void DrawBallf(bool * grid, uint16_t pwidth, uint16_t pheight, float x0, float y0, float z0, float radius) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 - r, t, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 + r, t, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 - t, r, radiusError);
		DrawDiskXY(grid, pwidth, pheight, x0, y0, z0 + t, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 - t, y0, z0, r, radiusError);
		DrawCircleYZ(grid, pwidth, pheight, x0 + t, y0, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 + t, z0, r, radiusError);
		DrawCircleXZ(grid, pwidth, pheight, x0, y0 - t, z0, r, radiusError);
//		DrawDiskYZ(grid, pwidth, pheight, x0 - t, y0, z0, r, radiusError);
//		DrawDiskYZ(grid, pwidth, pheight, x0 + t, y0, z0, r, radiusError);
//		DrawDiskYZ(grid, pwidth, pheight, x0 - r, y0, z0, t, radiusError);
//		DrawDiskYZ(grid, pwidth, pheight, x0 + r, y0, z0, t, radiusError);
//		DrawDiskXZ(grid, pwidth, pheight, x0, y0 + t, z0, r, radiusError);
//		DrawDiskXZ(grid, pwidth, pheight, x0, y0 - t, z0, r, radiusError);
//		DrawDiskXZ(grid, pwidth, pheight, x0, y0 + r, z0, t, radiusError);
//		DrawDiskXZ(grid, pwidth, pheight, x0, y0 - r, z0, t, radiusError);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			--r;
			radiusError += 2 * (t - r) + 1;
		}
	}
};

inline static void DrawBall(int radius, vector<voxel_offset> & ball_offsets) {
	int r = radius, t = 0;
	int radiusError = 1 - r;

	while (r >= t) {
		// pass in base point (x0,y0,z0), this algorithm's t as the current plane,
		// this algorithm's r as the radius, and pass along radius error.
		DrawDiskXY(-r, t, radiusError, ball_offsets);
		DrawDiskXY( r, t, radiusError, ball_offsets);
		DrawDiskXY(-t, r, radiusError, ball_offsets);
		DrawDiskXY( t, r, radiusError, ball_offsets);
		DrawCircleYZ(-t, r, radiusError, ball_offsets);
		DrawCircleYZ( t, r, radiusError, ball_offsets);
		DrawCircleXZ( t, r, radiusError, ball_offsets);
		DrawCircleXZ(-t, r, radiusError, ball_offsets);
		++t;
		if (radiusError < 0) {
			radiusError += 2 * t + 1;
		} else {
			r--;
			radiusError += 2 * (t - r) + 1;
		}
	}
};
#endif /* SPHERE_DRAWSPHERE_H_ */
