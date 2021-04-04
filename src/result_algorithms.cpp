/*
 * result_algorithms.cpp
 *
 *  Created on: 4 abr. 2021
 *      Author: Ruben Girela Castellón
 */


#include "../include/result_algorithms.h"

int ResultAlgorithms::Infeasable(vector<vector<int>> clusters, mat matriz){
	int restrictions = 0, col = 0, row=0;
	//walk through each cluster
	for(vector<vector<int>>::iterator it = clusters.begin(); it != clusters.end(); ++it){
		//walk through each node of the cluster
		for(vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
			//walk through nodes of matrix
			for(unsigned int i = 0; i<matriz.n_cols; ++i){
				//Is a triangular matrix
				if(i<(unsigned int)(*it2)){
					col = i;
					row = (*it2);
				}else{
					col = (*it2);
					row = i;
				}

				//if it isn't the identity
				if(col != row){

					//if the restriction is CL and the node is in the current cluster
					if(matriz(col,row)== -1 and find((*it).begin(),(*it).end(),i) != (*it).end())
							++restrictions;//increased constraint
					//in case it isn't and the constraint is ML
					else if(matriz(col,row) == 1 and find((*it).begin(),(*it).end(),i) == (*it).end())
							++restrictions;//increased constraint

				}
			}
		}
	}
	return restrictions;
}

