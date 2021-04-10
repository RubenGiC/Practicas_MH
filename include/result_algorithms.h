/*
 * result_algorithms.h
 *
 *  Created on: 4 abr. 2021
 *      Author: Ruben Girela Castellón
 */

#ifndef INCLUDE_RESULT_ALGORITHMS_H_
#define INCLUDE_RESULT_ALGORITHMS_H_

#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <armadillo>
using namespace std;
using namespace arma;

class ResultAlgorithms{

public:
	int Infeasable(vector<pair<int,int>> ML, vector<pair<int,int>> CL, vector<int> S);
	float Distance(vector<vector<int>> clusters, vector<vector<float>> atributos,vector<vector<float>> centroides);
	float distanciaEuclidea(vector<float> nod1, vector<float> nod2);
	float Fitness(vector<vector<float>> atributos, mat matriz, vector<vector<int>> clusters, vector<vector<float>> centroides, vector<pair<int,int>> ML, vector<pair<int,int>> CL, vector<int> S);
	//calculate Landa
	float createLanda(vector<vector<float>> atributos, mat matriz);
	//calculate the general deviation
	float generalDeviation(vector<vector<int>> v_clust, vector<vector<float>> atributos, vector<vector<float>> centroides);
};

#endif /* INCLUDE_RESULT_ALGORITHMS_H_ */
