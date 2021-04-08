/*
 * result_algorithms.cpp
 *
 *  Created on: 4 abr. 2021
 *      Author: Ruben Girela Castellón
 */


#include "../include/result_algorithms.h"

int ResultAlgorithms::Infeasable(vector<vector<int>> clusters, vector<pair<int,int>> ML, vector<pair<int,int>> CL){
	int restrictions = 0, col = 0, row=0;
	bool exit = false;

	for (unsigned int i = 0; i<ML.size(); ++i) {
		exit = false;
		for(unsigned int e = 0; e < clusters.size() && !exit; ++e){
			if(find(clusters[e].begin(), clusters[e].end(),ML[i].first) != clusters[e].end()){
				exit = true;
				if(find(clusters[e].begin(), clusters[e].end(),ML[i].second) == clusters[e].end()){
					++restrictions;
				}
			}else if(find(clusters[e].begin(), clusters[e].end(),ML[i].second) != clusters[e].end()){
				exit = true;
				if(find(clusters[e].begin(), clusters[e].end(),ML[i].first) == clusters[e].end()){
					++restrictions;
				}
			}
		}
	}

	for (unsigned int i = 0; i<CL.size(); ++i) {
			exit = false;
			for(unsigned int e = 0; e < clusters.size() && !exit; ++e){
				if(find(clusters[e].begin(), clusters[e].end(),CL[i].first) != clusters[e].end()){
					exit = true;
					if(find(clusters[e].begin(), clusters[e].end(),CL[i].second) != clusters[e].end()){
						++restrictions;
					}
				}else if(find(clusters[e].begin(), clusters[e].end(),CL[i].second) != clusters[e].end()){
					exit = true;
					if(find(clusters[e].begin(), clusters[e].end(),CL[i].first) != clusters[e].end()){
						++restrictions;
					}
				}
			}
		}

	/*for (const auto &it: ML) {
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
		}*/

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

//calculate the fitness
float ResultAlgorithms::Fitness(vector<vector<float>> atributos, mat matriz, vector<vector<int>> clusters, vector<vector<float>> centroides){

	vector<int> S(atributos.size());

	for(unsigned int e = 0; e < clusters.size(); ++e){
		for(vector<int>::iterator it = clusters[e].begin(); it != clusters[e].end(); ++it)
			S[*it] = e;
	}

	float f = generalDeviation(clusters, atributos, centroides) + (infeasibility(S, matriz) * createLanda(atributos, matriz));

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

	lan = lan/rest.size();

	return lan;
}

//calculate the general deviation
float ResultAlgorithms::generalDeviation(vector<vector<int>> v_clust, vector<vector<float>> atributos, vector<vector<float>> centroides){
	float distance = 0;

	//walk through each cluster
	for(unsigned int i = 0; i< v_clust.size(); ++i){
		//sumatorry(euclidean distance of all nodes in the cluster)
		for(vector<int>::iterator it = v_clust[i].begin(); it != v_clust[i].end(); ++it){
			distance += distanciaEuclidea(atributos[(*it)],centroides[i]);
		}
		//mean intra-cluster distance
		distance = distance / v_clust[i].size();
	}

	return distance/v_clust.size();
}

//calculate infeasibility
int ResultAlgorithms::infeasibility(vector<int> S, mat matriz){
	int restrictions = 0;

	//walk through each restriction
	for(unsigned int i = 1; i < matriz.n_cols; ++i)
		for(unsigned int e = 0; e < i; ++e)
			/*
			 * if the restriction is CL and the node is in the current cluster
			 * or
			 * it isn't and the constraint is ML
			*/
			if((matriz(i,e) > 0 && S[i] != S[e]) || (matriz(i,e) < 0 && S[i] == S[e]))
				++restrictions;

	return restrictions;
}
