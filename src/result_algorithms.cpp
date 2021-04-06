/*
 * result_algorithms.cpp
 *
 *  Created on: 4 abr. 2021
 *      Author: Ruben Girela Castellón
 */


#include "../include/result_algorithms.h"

int ResultAlgorithms::Infeasable(vector<vector<int>> clusters, vector<pair<int,int>> ML, vector<pair<int,int>> CL){
	int restrictions = 0, col = 0, row=0;

	for (const auto &it: ML) {
		for(const auto &it2: clusters){
			if(find(it2.begin(),it2.end(), it.first) != it2.end()){
				if(find(it2.begin(),it2.end(), it.second) == it2.end()){
					//cout << "Incumple ML: " << it.first << ", " << it.second << endl;
					++restrictions;
				}
			}
		}
	}

	for (const auto &it: CL) {
			for(const auto &it2: clusters){
				if(find(it2.begin(),it2.end(), it.first) != it2.end()){
					if(find(it2.begin(),it2.end(), it.second) != it2.end()){
						//cout << "Incumple CL: " << it.first << ", " << it.second << endl;
						++restrictions;
					}
				}
			}
		}

	/*//walk through each cluster
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
	}*/
	return restrictions;
}
//calculate the general deviation
float ResultAlgorithms::Distance(vector<vector<int>> clusters, vector<vector<float>> atributos,vector<vector<float>> centroides){
	float distance = 0;
	float intracluster = 0;
	//walk through each cluster
	for(unsigned int i = 0; i< clusters.size(); ++i){
		intracluster = 0;
		//sumatorry(euclidean distance of all nodes in the cluster)
		for(vector<int>::iterator it = clusters[i].begin(); it != clusters[i].end(); ++it){
			intracluster += distanciaEuclidea(atributos[(*it)],centroides[i]);
		}
		//mean intra-cluster distance
		intracluster = intracluster / clusters[i].size();
		distance += intracluster;
	}

	return distance/clusters.size();
}

//calcula la distancia euclidea entre 2 nodos
float ResultAlgorithms::distanciaEuclidea(vector<float> nod1, vector<float> nod2){
	//la formula es: sqrt(sumatoria((a_i - b_i)²))
	float suma=0;

	//sumatoria((a_i - b_i)²)
	for(unsigned int i=0; i<nod1.size();++i){
		suma += pow((nod1[i]-nod2[i]),2);
	}
	//sqrt(sumatoria)
	return sqrt(suma);
}
