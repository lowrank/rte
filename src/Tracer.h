/*
 * Tracer.h
 *
 * For MATLAB mex use. Caution: requires large memory.
 *
 *  Created on: Sep 14, 2017
 *      Author: lurker
 */

#ifndef Tracer_H
#define Tracer_H


#include <cstddef>
#include <iterator>
#include <memory>
#include <ostream>
#include <iostream>
#include <set>
#include <tuple>
#include <type_traits>
#include <utility>
#include <valarray>
#include <cstring>
#include <vector>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <set>
#include <queue>
#include <mexplus.h>
#include "common.h"

using std::vector;
using std::swap;
using std::cout;
using std::endl;

typedef struct Raylet {
	int32_t elem;
	double first[2];
	double second[2];
} Raylet;


class Tracer {
public:
	Tracer(M_Ptr angle, M_Ptr node, M_Ptr elem, M_Ptr neigh);
	virtual ~Tracer();


	void RayIntHelper(size_t numberofelems, size_t numberofnodesperelem,
			size_t numberofnodes,
			int32_t* pelems, double* pnodes, int32_t* pneighbors,
			int32_t i, double theta);

	void RayTrace(std::vector<double>& tmp, bool& intersect, double& q_t, double& q_eta,
			double& q_x1, double& q_y1, double& q_x2,
			double& q_y2, double& q_x3, double& q_y3,
			double& a, double& b, double& theta);

	void RayTrim(std::vector<double>& tmp, double &a, double &b);

	void RayShow();

	vector<vector<vector<Raylet>>> Ray;

};

inline double INTERSECT_DET(double X1, double Y1, double X2, double Y2, double THETA){
		return ((cos((THETA)) * ((Y1) - (Y2))) - (sin((THETA)) * ((X1) - (X2))));
	}
inline double INTERSECT_CROSS(double X1,double Y1,double X2,double Y2,double A,double B) {
		return ((((X1) - (X2)) * ((B) - (Y2))) -(((A) - (X2)) * ((Y1)- (Y2))));
	}


#endif /* Tracer.h_H */
