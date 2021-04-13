/*
 * result_algorithms.cpp
 *
 *  Created on: 4 abr. 2021
 *      Author: Ruben Girela Castellón
 */


#include "../include/result_algorithms.h"
#include <string>
#include <stdlib.h> //abs


//calculates the number of constrains it violates
int ResultAlgorithms::Infeasable(vector<pair<int,int>> ML, vector<pair<int,int>> CL, vector<int> S){
	int restrictions = 0;
	unsigned int max = CL.size();;

	//get the maximum size
	if(ML.size() > CL.size())
		max = ML.size();

	//walk through each constrain ML and CL
	for (unsigned int i = 0; i<max; ++i) {

		//if i is less than the size of ML
		if(i < ML.size()){
			//if the nodes have diferent clusters
			if(S[ML[i].first] != S[ML[i].second]){
				//violates the constrain
				++restrictions;
			}
		}
		//if i is less than the size of CL
		if(i < CL.size()){
			//if the nodes have the same clusters
			if(S[CL[i].first] == S[CL[i].second]){
				//violates the constrain
				++restrictions;
			}
		}
	}

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

//calculate the distance error
float ResultAlgorithms::ErrorDistance(vector<vector<int>> clusters, vector<vector<float>> atributos,vector<vector<float>> centroides, string type_data_file){
	float error_distance = 0, optimal_distance = 0;
	//calculate the distance
	error_distance = Distance(clusters, atributos, centroides);

	//choose the optimal distance
	if(type_data_file.compare("ZOO") == 0)
		optimal_distance = 0.904799856;
	else if(type_data_file.compare("GLASS") == 0)
		optimal_distance = 0.364290282;
	else
		optimal_distance = 0.220423749;

	//calculate the distance error
	error_distance = abs(error_distance-optimal_distance);

	return error_distance;
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

//calculate the fitness
float ResultAlgorithms::Fitness(vector<vector<float>> atributos, mat matriz, vector<vector<int>> clusters, vector<vector<float>> centroides, vector<pair<int,int>> ML, vector<pair<int,int>> CL, vector<int> S){
	//fitness = general deviation + (infeasable * landa)
	float f = generalDeviation(clusters, atributos, centroides) + (Infeasable(ML, CL, S) * createLanda(atributos, matriz));

	return f;
}

//calculate Landa
float ResultAlgorithms::createLanda(vector<vector<float>> atributos, mat matriz){
	float lan = 0, actual_distance=0;

	//calculate the maximum distance
	for(unsigned int i = 0; i < atributos.size(); ++i){
		for(unsigned int e = i+1; e < atributos.size(); ++e){
			actual_distance = distanciaEuclidea(atributos[i],atributos[e]);
			if(actual_distance > lan){
				lan = actual_distance;
			}
		}
	}
	//count the number of restrictions total
	uvec rest = (find(matriz == 1 or matriz == -1));

	//max distance / total number of problem restrictions
	lan = lan/rest.size();

	return lan;
}

//calculate the general deviation
float ResultAlgorithms::generalDeviation(vector<vector<int>> v_clust, vector<vector<float>> atributos, vector<vector<float>> centroides){
	float distance = 0, intra_cluster = 0;

	//walk through each cluster
	for(unsigned int i = 0; i< v_clust.size(); ++i){
		distance = 0;
		//sumatorry(euclidean distance of all nodes in the cluster)
		for(vector<int>::iterator it = v_clust[i].begin(); it != v_clust[i].end(); ++it){
			distance += distanciaEuclidea(atributos[(*it)],centroides[i]);
		}
		//mean intra-cluster distance
		intra_cluster += distance / v_clust[i].size();
	}

	return intra_cluster/v_clust.size();
}
