/*
 * Tracer.cpp
 *
 *  Created on: Sep 14, 2017
 *      Author: lurker
 */

#include "Tracer.h"

Tracer::Tracer(M_Ptr angle, M_Ptr nodes, M_Ptr elems, M_Ptr neigh){
	auto numberofangle = mxGetN(angle);
	auto numberofnodes = mxGetN(nodes);
	auto numberofelems = mxGetN(elems);
	auto numberofnodesperelem = mxGetM(elems);

	auto pangle = mxGetPr(angle);
	auto pnodes = mxGetPr(nodes);
	auto pelems = (int*)mxGetPr(elems);
	auto pneigh = (int*)mxGetPr(neigh);


	Ray.resize(numberofangle);
#pragma omp parallel for
	for (int aId = 0; aId < numberofangle; ++aId) {
		Ray[aId].resize(numberofnodes);
		double theta = pangle[aId];

		RayIntHelper(numberofelems, numberofnodesperelem, numberofnodes, pelems, pnodes, pneigh, aId, theta);
	}
}
Tracer::~Tracer() {}

void Tracer::RayIntHelper(size_t numberofelems, size_t numberofnodesperelem,
		size_t numberofnodes,
		int32_t* pelems, double* pnodes, int32_t* pneighbors,int32_t i, double theta){


	mwSize vertex, e_vl, e_vr;
	double x1, x2, y1, y2,  a, b;
	double ret, t, eta;
	double q_x1, q_y1, q_x2, q_y2, q_x3, q_y3;
	double q_eta, q_t;
	int32_t index;
	bool intersect;
	Raylet tmp_ray;

	std::vector<double> tmp;
	std::unordered_set<int> s_elem;
	std::queue<int> q_elem;
	std::unordered_set<int> s_node;

	/*
	 * for loop, main part
	 */
	for (int j = 0; j < numberofelems; j++) {
		for (int k = 0; k < numberofnodesperelem; k++){
			/*
			 * index of vertex in nodes.
			 */
			vertex = pelems[numberofnodesperelem * j + k] - 1;

			if (s_node.find(vertex) != s_node.end()) {
				continue;
			}

			s_node.insert(vertex);
			/*
			 * calculate the ray intersects with the boundary.
			 */
			a = pnodes[2 * vertex    ];
			b = pnodes[2 * vertex + 1];

			// clean up set.
			s_elem.clear();
			// push current element
			q_elem.push(j);
			s_elem.insert(j);
			/*
			 * check the jth.
			 */
			index = j;

			q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
			q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

			q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
			q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

			q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
			q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];

			RayTrace(tmp, intersect, q_t, q_eta,
					q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
					a, b, theta);

			RayTrim(tmp, a, b);

			if (tmp.size()) {
					tmp_ray.elem = index;
				if (tmp.size() == 2) {
					tmp_ray.first[0] = a;
					tmp_ray.first[1] = b;
					tmp_ray.second[0] = tmp[0];
					tmp_ray.second[1] = tmp[1];
				}
				else {
					tmp_ray.first[0] = tmp[0];
					tmp_ray.first[1] = tmp[1];
					tmp_ray.second[0] = tmp[2];
					tmp_ray.second[1] = tmp[3];
				}

				/*
				 * remove duplicates, since sorted, only test first two.
				 */

				if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
					fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){

					Ray[i][vertex].push_back(tmp_ray);
				}
			}
			while (!q_elem.empty()){
				auto top = q_elem.front();
				q_elem.pop();
				// loop over neighbors.
				// fail safe strategy, if it is possible to intersect, must push.
				for (int32_t q_elem_i = 0; q_elem_i < 3; q_elem_i ++){

					index = pneighbors[3 * (top) + q_elem_i] - 1;

					// test if it is going to push them in.
					// first : valid index of element
					if (index > -1 && s_elem.find(index) == s_elem.end()) {
						/*
						 * calculate the open angle limit. And find the suitable one.
						 *
						 * Target function is max(min(abs(langle - angle), abs(rangle - angle)))
						 *
						 * if both gives 0, then push both into queue.
						 */

						q_x1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1)];
						q_y1 = pnodes[2 * (pelems[numberofnodesperelem * index] - 1) + 1];

						q_x2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1)];
						q_y2 = pnodes[2 * (pelems[numberofnodesperelem * index + 1] - 1) + 1];

						q_x3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1)];
						q_y3 = pnodes[2 * (pelems[numberofnodesperelem * index + 2] - 1) + 1];


						/*
						 * [a , b] 's angle for the element,
						 */

						RayTrace(tmp, intersect, q_t, q_eta,
								q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
								a, b, theta);

						if (intersect){
							q_elem.push(index);
							s_elem.insert(index);

							/*
							 * calculate integral over this element
							 *
							 * first step: calculate the intersection,
							 * second step: integral
							 */
							RayTrim(tmp, a, b);

							if (tmp.size()) {
									tmp_ray.elem = index;
								if (tmp.size() == 2) {
									tmp_ray.first[0] = a;
									tmp_ray.first[1] = b;
									tmp_ray.second[0] = tmp[0];
									tmp_ray.second[1] = tmp[1];
								}
								else {
									tmp_ray.first[0] = tmp[0];
									tmp_ray.first[1] = tmp[1];
									tmp_ray.second[0] = tmp[2];
									tmp_ray.second[1] = tmp[3];
								}

								/*
								 * remove duplicates, since sorted, only test first two.
								 */
								if (!Ray[i][vertex].size() || fabs(tmp_ray.first[0] - Ray[i][vertex].back().first[0]) > MEX_EPS ||
									fabs(tmp_ray.first[1] - Ray[i][vertex].back().first[1]) > MEX_EPS){
									Ray[i][vertex].push_back(tmp_ray);

								}
							}
						}
					}
				}
			}
			Ray[i][vertex].shrink_to_fit();
		}
	}
}

void Tracer::RayTrace(std::vector<double>& tmp, bool& intersect, double& q_t, double& q_eta,
			double& q_x1, double& q_y1, double& q_x2,
			double& q_y2, double& q_x3, double& q_y3,
			double& a, double& b, double& theta){

	intersect = false;
	tmp.clear();
	// 1 - 2 intersects.
	if (fabs(INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta)) < 1 * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)) < 1 * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x1 - q_x2) + MEX_EPS < fabs(q_x1 - a) + fabs(q_x2 - a) ||
					fabs(q_y1 - q_y2) + MEX_EPS < fabs(q_y1 - b) + fabs(q_y2 - b) ){
				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x1, q_y1, q_x2, q_y2, a, b)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);

		q_eta = INTERSECT_DET(a, b, q_x2, q_y2, theta)/
				INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x1 + (1- q_eta) * q_x2);
			tmp.push_back(q_eta * q_y1 + (1 -q_eta) * q_y2);
		}
	}
	// 2 - 3
	if (fabs(INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta)) < 1 * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)) < 1 * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x2 - q_x3) + MEX_EPS < fabs(q_x2 - a) + fabs(q_x3 - a) ||
					fabs(q_y2 - q_y3) + MEX_EPS < fabs(q_y2 - b) + fabs(q_y3 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x2, q_y2, q_x3, q_y3, a, b)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		q_eta = INTERSECT_DET(a, b, q_x3, q_y3, theta)/
				INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta);

		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x2 + (1 - q_eta) * q_x3);
			tmp.push_back(q_eta * q_y2 + (1 - q_eta) * q_y3);
		}
	}
	// 3 - 1
	if (fabs(INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta)) < 1 * MEX_EPS){
		// ray parallel for edge.
		if (fabs(INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)) < 1 * MEX_EPS){
			// it is lucky to be colinear.
			if (fabs(q_x3 - q_x1) + MEX_EPS < fabs(q_x3 - a) + fabs(q_x1 - a) ||
					fabs(q_y3 - q_y1) + MEX_EPS < fabs(q_y3 - b) + fabs(q_y1 - b) ){

				// outside
			}
			else {
				intersect = true;
				tmp.push_back(a);
				tmp.push_back(b);
			}
		}
		else {
			// intersect = false;
		}
	}
	else{
		// not parallel, then there is a intersect.
		q_t =INTERSECT_CROSS(q_x3, q_y3, q_x1, q_y1, a, b)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);

		q_eta = INTERSECT_DET(a, b, q_x1, q_y1, theta)/
				INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta);
		/*
		 * q_t = 0 means colinear
		 *
		 */
		if (q_t >= 0 && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
			intersect = true;
			tmp.push_back(q_eta * q_x3 + (1 - q_eta) * q_x1);
			tmp.push_back(q_eta * q_y3 + (1 - q_eta) * q_y1);
		}
	}
}

void Tracer::RayTrim(std::vector<double>& tmp, double &a, double &b){
	for (int32_t tmp_i = 0; tmp_i < tmp.size()/2; tmp_i ++){

		if (fabs(a - tmp[2 * tmp_i]) + fabs(b - tmp[2 * tmp_i + 1]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2 * tmp_i, tmp.begin() + 2 * tmp_i + 2);
			tmp_i --;
		}
	}

	if (tmp.size() == 6) {
	// 2 - 3 duplicates
		if (fabs(tmp[2] - tmp[4]) + fabs(tmp[3] - tmp[5]) < MEX_EPS){
			tmp.erase(tmp.begin() + 2, tmp.begin() + 4);
		}
	// 3 - 1 duplicates or 1 - 2 duplicates
		else {
			tmp.erase(tmp.begin(), tmp.begin() + 2);
		}
	}
	if (tmp.size() == 4) {
		//remove duplicates
		if (fabs(tmp[0] - tmp[2]) + fabs(tmp[1] - tmp[3]) < MEX_EPS) {
			tmp.clear();
		}
		else {
			// reorder
			if (pow(tmp[0] - a, 2) + pow(tmp[1] - b, 2) > pow(tmp[2] - a, 2) + pow(tmp[3] - b, 2)){
				swap(tmp[0], tmp[2]);
				swap(tmp[1], tmp[3]);
			}
		}
	}
}

void Tracer::RayShow(){
	int32_t tmp_i, tmp_j;
	size_t tmp_total = 0;
	if (Ray.size() != 0) {
		for (int32_t i = 0 ; i < Ray.size(); i++){
			tmp_i = Ray[i].size();
			for (int32_t j = 0; j < tmp_i; j++){
				tmp_j = Ray[i][j].capacity();
				tmp_total += tmp_j * 40;

//				for (int32_t k = 0; k < tmp_j; k++) {
//					std::cout << i << "th Angle, "
//							<< j << "th node, "
//							<< k << "th raylet: passes through "
//							<< Ray[i][j][k].elem << ", starting from "
//							<< Ray[i][j][k].first[0] << ", " << Ray[i][j][k].first[1] << " --> "
//							<< Ray[i][j][k].second[0] << ", " << Ray[i][j][k].second[1] << std::endl;
//				}
			}
		}
	}
	std::cout
	<< tmp_total / 1024.0/ 1024.0/ 1024.0
	<< " GBytes used in geometry storage."
	<< std::endl;
}

using namespace mexplus;

template class Session<Tracer>;

namespace {
    MEX_DEFINE(new)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 4);
        OutputArguments output(nlhs, plhs, 1);

        output.set(0, Session<Tracer>::create(new Tracer(C_CAST(prhs[0]), C_CAST(prhs[1]), C_CAST(prhs[2]), C_CAST(prhs[3]))));

    }

    MEX_DEFINE(delete)(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        InputArguments input(nrhs, prhs, 1);
        OutputArguments output(nlhs, plhs, 0);
        Session<Tracer>::destroy(input.get(0));
    }

    MEX_DEFINE(disp) (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    	InputArguments input(nrhs, prhs, 1);
    	OutputArguments output(nlhs, plhs, 0);

    	auto DOM = Session<Tracer>::get(input.get(0));

    	DOM->RayShow();
    }



}

MEX_DISPATCH


