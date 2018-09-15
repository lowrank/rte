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

        RayHelper(numberofelems, numberofnodesperelem, numberofnodes, pelems, pnodes, pneigh, aId, theta);
    }
}
Tracer::~Tracer() {}

void Tracer::RayHelper(size_t numberofelems, size_t numberofnodesperelem,
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
     *
     *
     * s_node stores the vertex visited for this ray.
     * s_elem stores the element visited.
     *
     * a, b is the starting point for current ray.
     * theta as the direction angle.
     *
     * q_elem is the temp queue for testing neighbors.
     *
     *
     * tmp is the output for this run.
     *
     */
    for (int j = 0; j < numberofelems; j++) {
        for (int k = 0; k < numberofnodesperelem; k++){



            vertex = pelems[numberofnodesperelem * j + k] - 1;  // index of vertex in nodes of jth element. minus 1 to be 0-based.

            if (s_node.find(vertex) != s_node.end()) {
                continue;                                       // make sure that we have not run for this vertex yet.
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
             * check the intersection inside jth element.
             */

            q_x1 = pnodes[2 * (pelems[numberofnodesperelem * j] - 1)];
            q_y1 = pnodes[2 * (pelems[numberofnodesperelem * j] - 1) + 1];

            q_x2 = pnodes[2 * (pelems[numberofnodesperelem * j + 1] - 1)];
            q_y2 = pnodes[2 * (pelems[numberofnodesperelem * j + 1] - 1) + 1];

            q_x3 = pnodes[2 * (pelems[numberofnodesperelem * j + 2] - 1)];
            q_y3 = pnodes[2 * (pelems[numberofnodesperelem * j + 2] - 1) + 1];


            /*
             * find possible intersection with the ray and the triangle.
             * it is possible to have 0, 1, 2 intersections depending on various cases.
             * Just compute the intersection of the ray with 3 sides.
             * We got 3 intersections. We arrange them according to the distance from (a,b) and abandon the duplicates.
             *
             */
            RayTrace(tmp, intersect, q_t, q_eta,
                     q_x1, q_y1, q_x2, q_y2, q_x3, q_y3,
                     a, b, theta);

            RayTrim(tmp, a, b);

            if (tmp.size()) {
                tmp_ray.elem = j;
                if (tmp.size() == 4) {
                    tmp_ray.first[0] = tmp[0];
                    tmp_ray.first[1] = tmp[1];
                    tmp_ray.second[0] = tmp[2];
                    tmp_ray.second[1] = tmp[3];
                }
                else if (tmp.size() == 2) {
                    tmp_ray.first[0] = a;
                    tmp_ray.first[1] = b;
                    tmp_ray.second[0] = tmp[0];
                    tmp_ray.second[1] = tmp[1];
                }

                Ray[i][vertex].push_back(tmp_ray);

            }

            // if queue is not empty, means we have not tried all elements on this ray
            while (!q_elem.empty()){
                auto top = q_elem.front();
                q_elem.pop();
                // loop over neighbors.

                for (auto nnid = 0; nnid < 3; ++nnid){

                    index = pneighbors[3 * top + nnid] - 1;

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

                                if (tmp.size() == 4) {
                                    tmp_ray.first[0] = tmp[0];
                                    tmp_ray.first[1] = tmp[1];
                                    tmp_ray.second[0] = tmp[2];
                                    tmp_ray.second[1] = tmp[3];
                                }
                                else if (tmp.size() == 2) {
                                    tmp_ray.first[0] = a;
                                    tmp_ray.first[1] = b;
                                    tmp_ray.second[0] = tmp[0];
                                    tmp_ray.second[1] = tmp[1];
                                }

                                /*
                                 * remove duplicates, since sorted, only test first two.
                                 */
                                Ray[i][vertex].push_back(tmp_ray);

                            }
                        }
                    }
                }
            }

            auto cmp = [&](Raylet& A, Raylet& B) ->bool {
                return pow((A.first[0] + A.second[0])/2 - a, 2) + pow((A.first[1] + A.second[1])/2 - b, 2)
                       < pow((B.first[0] + B.second[0])/2 - a, 2) + pow((B.first[1] + B.second[1])/2 - b, 2);
            };

            sort(Ray[i][vertex].begin(), Ray[i][vertex].end(), cmp);


            for (auto it = Ray[i][vertex].begin(); it != Ray[i][vertex].end(); ) {

                if ( fabs(it->first[0] - it->second[0]) < MEX_EPS
                     && fabs(it->first[1] - it->second[1]) < MEX_EPS ) {
                    it = Ray[i][vertex].erase(it);
                }
                else {
                    ++it;
                }

            }

            auto itrFirst = Ray[i][vertex].begin();
            auto itrLast = Ray[i][vertex].end();
            auto itrResult = itrFirst;
            if (itrFirst == itrLast) { // only one raylet.
                itrResult = itrLast;
            }
            else {
                itrResult = itrFirst;
                while (++itrFirst != itrLast) {
                    if (!(fabs(itrResult->first[0] - itrFirst->first[0]) < MEX_EPS
                          && fabs(itrResult->first[1] - itrFirst->first[1]) < MEX_EPS)) {
                        // only test equal condiion for the first coordinate.
                        ++itrResult;
                        (itrResult)->elem = itrFirst->elem;
                        (itrResult)->first[0] = itrFirst->first[0];
                        (itrResult)->second[0] = itrFirst->second[0];
                        (itrResult)->first[1] = itrFirst->first[1];
                        (itrResult)->second[1] = itrFirst->second[1];
                    }
                }
                ++itrResult;
            }
            Ray[i][vertex].resize( std::distance(Ray[i][vertex].begin(),itrResult) );
            Ray[i][vertex].shrink_to_fit();
        }
    }

}

/*
 *
 *  RayTrace finds the intersection point for a ray and 3 segments.
 *
 *  there will be 3 intersections at most (with duplicates). The triangle is not degenerated.
 *
 *
 */
void Tracer::RayTrace(std::vector<double>& tmp, bool& intersect, double& q_t, double& q_eta,
                      double& q_x1, double& q_y1, double& q_x2,
                      double& q_y2, double& q_x3, double& q_y3,
                      double& a, double& b, double& theta){

    intersect = false;
    tmp.clear();
    // 1 - 2 intersects.
    if (fabs(INTERSECT_DET(q_x1, q_y1, q_x2, q_y2, theta)) < 1 * MEX_EPS){
        // ray parallel for edge.
        // in this case, we consider the intersection is not valid.
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
        if (q_t >= -MEX_EPS && q_eta >=-MEX_EPS && q_eta <= 1 + MEX_EPS) {
            intersect = true;
            tmp.push_back(q_eta * q_x1 + (1- q_eta) * q_x2);
            tmp.push_back(q_eta * q_y1 + (1 -q_eta) * q_y2);
        }
    }
    // 2 - 3
    if (fabs(INTERSECT_DET(q_x2, q_y2, q_x3, q_y3, theta)) < 1 * MEX_EPS){
        // ray parallel for edge.
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
        if (q_t >= -MEX_EPS && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
            intersect = true;
            tmp.push_back(q_eta * q_x2 + (1 - q_eta) * q_x3);
            tmp.push_back(q_eta * q_y2 + (1 - q_eta) * q_y3);
        }
    }
    // 3 - 1
    if (fabs(INTERSECT_DET(q_x3, q_y3, q_x1, q_y1, theta)) < 1 * MEX_EPS){
        // ray parallel for edge.
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
        if (q_t >= MEX_EPS && q_eta >= -MEX_EPS && q_eta <= 1 + MEX_EPS) {
            intersect = true;
            tmp.push_back(q_eta * q_x3 + (1 - q_eta) * q_x1);
            tmp.push_back(q_eta * q_y3 + (1 - q_eta) * q_y1);
        }
    }
    // the output has at most 3 points.
}

void Tracer::RayTrim(std::vector<double>& tmp, double &a, double &b){

    assert(tmp.size() <= 6);
    if (tmp.size() == 0) return;
    /*
     * there are 3 points at most, to reduce the number, we have to find out which to be reduced.
     *
     * for the case with 3 points (size = 6), there is a duplicate.
     *
     */
    if (tmp.size() == 6) {
        // 2 - 3 duplicates
        if (fabs(tmp[2] - tmp[4]) + fabs(tmp[3] - tmp[5]) < MEX_EPS){
            tmp.erase(tmp.begin() + 2, tmp.begin() + 4);
        }
            // 3 - 1 duplicates or 1 - 2 duplicates
        else if (fabs(tmp[0] - tmp[2]) + fabs(tmp[1] - tmp[3]) < MEX_EPS) {
            tmp.erase(tmp.begin(), tmp.begin() + 2);
        }
        else if (fabs(tmp[0] - tmp[4]) + fabs(tmp[1] - tmp[5]) < MEX_EPS) {
            tmp.erase(tmp.begin(), tmp.begin() + 2);
        }
    }

    if (tmp.size() == 4) {
        // check duplicates

        // reorder
        if (pow(tmp[0] - a, 2) + pow(tmp[1] - b, 2) > pow(tmp[2] - a, 2) + pow(tmp[3] - b, 2)){
            swap(tmp[0], tmp[2]);
            swap(tmp[1], tmp[3]);
        }

    }
}
//TODO: get option for display verbose information.
void Tracer::RayShow(int l){
    if (l < 0) {
        int32_t tmp_i, tmp_j;
        size_t tmp_total = 0;
        if (Ray.size() != 0) {
            for (int32_t i = 0; i < Ray.size(); i++) {
                tmp_i = Ray[i].size();
                for (int32_t j = 0; j < tmp_i; j++) {
                    tmp_j = Ray[i][j].capacity();
                    tmp_total += tmp_j * sizeof(Raylet);

//                     for (int32_t k = 0; k < tmp_j; k++) {
// 
//                         std::cout << i << "th Angle, "
//                                   << j << "th node, "
//                                   << k << "th raylet: passes through "
//                                   << Ray[i][j][k].elem << ", starting from "
//                                   << Ray[i][j][k].first[0] << ", " << Ray[i][j][k].first[1] << " --> "
//                                   << Ray[i][j][k].second[0] << ", " << Ray[i][j][k].second[1] << std::endl;
// 
// 
//                     }
                }
            }
        }
        std::cout
                << tmp_total / 1024.0 / 1024.0 / 1024.0
                << " GBytes used in geometry storage."
                << std::endl;
    }
    else {
        int32_t tmp_i, tmp_j;
        size_t tmp_total = 0;
        if (Ray.size() != 0) {
            for (int32_t i = 0; i < Ray.size(); i++) {
                tmp_i = Ray[i].size();
                for (int32_t j = 0; j < tmp_i; j++) {
                    tmp_j = Ray[i][j].capacity();
                    tmp_total += tmp_j * sizeof(Raylet);

                    for (int32_t k = 0; k < tmp_j; k++) {
                        if (j == l) {
                            std::cout << i << "th Angle, "
                                      << j << "th node, "
                                      << k << "th raylet: passes through "
                                      << Ray[i][j][k].elem << ", starting from "
                                      << Ray[i][j][k].first[0] << ", " << Ray[i][j][k].first[1] << " --> "
                                      << Ray[i][j][k].second[0] << ", " << Ray[i][j][k].second[1] << std::endl;


                        }
                    }
                }
            }
        }
        std::cout
                << tmp_total / 1024.0 / 1024.0 / 1024.0
                << " GBytes used in geometry storage."
                << std::endl;
    }
}

void Tracer::RayBC(double* pnodes, size_t numberofnodes, int* pelems, size_t numberofnodesperelem, double* pval, double* sol, double* ptr){
    auto angles = Ray.size();
#pragma omp parallel for
    for (int si = 0; si < angles; ++si) { // angles
        double det, length, eta, lambda, lv, rv, acc, bd;
        double x1, y1, x2, y2, x3, y3;
        int vertex_1, vertex_2, vertex_3;
        for (int sj = 0; sj < numberofnodes; ++sj){ // nodes
            acc = 0.;
            if (Ray[si][sj].size()) { // if this coming direction is reversible
                /*
                 * get last element?
                 */
                for (auto it : Ray[si][sj]) { // traverse all raylets.
                    vertex_1 = pelems[it.elem * numberofnodesperelem] - 1;
                    vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
                    vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

                    x1 = pnodes[2 * vertex_1];
                    y1 = pnodes[2 * vertex_1 + 1];
                    x2 = pnodes[2 * vertex_2];
                    y2 = pnodes[2 * vertex_2 + 1];
                    x3 = pnodes[2 * vertex_3];
                    y3 = pnodes[2 * vertex_3 + 1];

                    det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
                    length = sqrt(pow(it.first[0] - it.second[0], 2) +
                                  pow(it.first[1] - it.second[1], 2));

                    eta = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
                    eta /= det;

                    lambda = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
                    lambda /= det;


                    lv = lambda * pval[vertex_1] + eta * pval[vertex_2] +
                         (1 - lambda - eta) * pval[vertex_3];

                    eta = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
                    eta /= det;

                    lambda = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
                    lambda /= det;


                    rv = lambda * pval[vertex_1] + eta * pval[vertex_2] +
                         (1 - lambda - eta) * pval[vertex_3];

                    acc += (rv + lv) * length / 2.0;
                } // all raylet

                bd = lambda * sol[angles * vertex_1 + si] + eta * sol[angles *  vertex_2 + si] +
                     (1 - lambda - eta) * sol[angles * vertex_3 + si];
                ptr[angles* sj + si] = exp(-acc) * bd;
            }// if trace
            else {
                ptr[angles * sj + si] = sol[angles * sj + si]; // only boundary term is used.
            }
        } // node
    } // angle
}


void Tracer::RayIN(double* pnodes, size_t numberofnodes, int* pelems, size_t numberofnodesperelem, double* pval, double* sol, double* ptr) {
    auto angles = Ray.size();
    for (int si = 0; si < angles; ++si) { // angles
        double det, length, eta1, lambda1, eta2, lambda2, lv, rv, acc, bd1, bd2;
        double x1, y1, x2, y2, x3, y3;
        int vertex_1, vertex_2, vertex_3;
        for (int sj = 0; sj < numberofnodes; ++sj) { // nodes
            acc = 0.;

            if (Ray[si][sj].size()) {
                for (auto it : Ray[si][sj]) {

                    vertex_1 = pelems[it.elem * numberofnodesperelem] - 1;
                    vertex_2 = pelems[it.elem * numberofnodesperelem + 1] - 1;
                    vertex_3 = pelems[it.elem * numberofnodesperelem + 2] - 1;

                    x1 = pnodes[2 * vertex_1];
                    y1 = pnodes[2 * vertex_1 + 1];
                    x2 = pnodes[2 * vertex_2];
                    y2 = pnodes[2 * vertex_2 + 1];
                    x3 = pnodes[2 * vertex_3];
                    y3 = pnodes[2 * vertex_3 + 1];

                    det = (x1 - x3) * (y2 - y3) - (x2 - x3) * (y1 - y3);
                    length = sqrt(pow(it.first[0] - it.second[0], 2) +
                                  pow(it.first[1] - it.second[1], 2));

                    // first node

                    eta1 = ((y3 - y1) * (it.first[0] - x3) + (x1 - x3) * (it.first[1] - y3));
                    eta1 /= det;

                    lambda1 = (y2 - y3) * (it.first[0] - x3) + (x3 -  x2) * (it.first[1] - y3);
                    lambda1 /= det;

                    // get values
                    lv = lambda1 * pval[vertex_1] + eta1 * pval[vertex_2] +
                         (1 - lambda1 - eta1) * pval[vertex_3];
                    bd1 = lambda1 * sol[angles * vertex_1 + si] + eta1 * sol[angles * vertex_2 + si] +
                         (1 - lambda1 - eta1) * sol[angles * vertex_3 + si];

                    // second node

                    eta2 = ((y3 - y1) * (it.second[0] - x3) + (x1 - x3) * (it.second[1] - y3));
                    eta2 /= det;

                    lambda2 = (y2 - y3) * (it.second[0] - x3) + (x3 -  x2) * (it.second[1] - y3);
                    lambda2 /= det;

                    // get values
                    rv = lambda2 * pval[vertex_1] + eta2 * pval[vertex_2] +
                         (1 - lambda2 - eta2) * pval[vertex_3];
                    bd2 = lambda2 * sol[angles * vertex_1 + si] + eta2 * sol[angles *  vertex_2 + si] +
                         (1 - lambda2 - eta2) * sol[angles * vertex_3 + si];

                    // inserting. Runge Kutta 3 stages
                    // 1st
                    ptr[angles* sj + si] += exp(-acc) * bd1 * length / 6.0;

                    acc += (0.5 * rv + 1.5 * lv) * length / 4.0;
                    // 2nd
                    ptr[angles* sj + si] += exp(-acc) * (bd1 + bd2) * length/ 3.0;

                    acc += (1.5 * rv + 0.5 * lv) * length / 4.0;
                    // 3rd
                    ptr[angles* sj + si] += exp(-acc) * bd2 * length/ 6.0;

                }
            }
        }// nodes
    }// angles

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
    InputArguments input(nrhs, prhs, 2);
    OutputArguments output(nlhs, plhs, 0);

    auto tracer = Session<Tracer>::get(input.get(0));
    auto pl      = mxGetPr(prhs[1]);

    tracer->RayShow(int(*pl));
}

MEX_DEFINE(boundary_transport)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    InputArguments input(nrhs, prhs, 5);
    OutputArguments output(nlhs, plhs, 1);

    auto tracer = Session<Tracer>::get(input.get(0));
    auto pnodes = mxGetPr(prhs[1]);
    auto numberofnodes = mxGetN(prhs[1]);

    auto pelems = (int*)mxGetPr(prhs[2]);
    auto numberofnodesperelem = mxGetM(prhs[2]);

    auto pval = mxGetPr(prhs[3]);
    auto sol = mxGetPr(prhs[4]);

    auto angles = tracer->Ray.size();

    plhs[0] = mxCreateNumericMatrix(angles, numberofnodes, mxDOUBLE_CLASS, mxREAL);
    auto ptr = mxGetPr(plhs[0]);

    tracer->RayBC(pnodes, numberofnodes,pelems, numberofnodesperelem, pval, sol, ptr);
}

MEX_DEFINE(interior_transport)(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
    InputArguments input(nrhs, prhs, 5);
    OutputArguments output(nlhs, plhs, 1);

    auto tracer = Session<Tracer>::get(input.get(0));
    auto pnodes = mxGetPr(prhs[1]);
    auto numberofnodes = mxGetN(prhs[1]);

    auto pelems = (int*)mxGetPr(prhs[2]);
    auto numberofnodesperelem = mxGetM(prhs[2]);

    auto pval = mxGetPr(prhs[3]);
    auto sol = mxGetPr(prhs[4]);

    auto angles = tracer->Ray.size();

    plhs[0] = mxCreateNumericMatrix(angles, numberofnodes, mxDOUBLE_CLASS, mxREAL);
    auto ptr = mxGetPr(plhs[0]);

    tracer->RayIN(pnodes, numberofnodes,pelems, numberofnodesperelem, pval, sol, ptr);
}


}

MEX_DISPATCH
