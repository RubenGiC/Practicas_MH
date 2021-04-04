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
	int Infeasable(vector<vector<int>> clusters, mat matriz);
};

#endif /* INCLUDE_RESULT_ALGORITHMS_H_ */
