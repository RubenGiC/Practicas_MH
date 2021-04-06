/*
 * lector_ficheros.cpp
 *
 *  Created on: 14 mar. 2021
 *      Author: Ruben Girela Castellón
 */
#include <iostream>
#include <armadillo>// para la matriz
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <math.h> // para calculos matematicos (sqrt, pow (x^y), ...)
#include <algorithm>// para barajar el vector
#include <ctime>// libreria de time para el shufle aleatorio de vectores
#include <stack>// use stack

#include "../include/PAR.h"
#include "../include/random.h"
using namespace std;
using namespace arma;

PAR::PAR(string fichero_set, string fichero_set_const){

	//creo el vector de atributos y la matriz de restricciones
	lectura(fichero_set, fichero_set_const);

	//en función del numero de atributos que tenga el numero de clusters sera 16 o 7
	switch(atributos[0].size()){
	case 5:
		k=16;
		break;
	default:
		k=7;
		break;
	}

	centroides.resize(k);//reservo el tamaño del vector de tamaño k
	clusters.resize(k);//reservo memoria para el vector de elementos de cada cluster

	Set_random(37);//creo una semilla para los valores aleatorios

	//Y genero aleatoriamente los centroides inicialmente distintos de dimensión n
	for(int i=0; i<k; ++i){
		for(unsigned int e=0; e <atributos[0].size(); ++e)
			centroides[i].push_back(Rand());
	}

	//insert all indices of atributos
	for(unsigned int i=0; i<atributos.size(); ++i){
		RSI.push_back(i);
	}

	//despues barajo los indices de los atributos
	srand(unsigned (37));//genero una semilla fija
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR::shuffleRSI(){
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR::lectura(string fichero_set, string fichero_set_const){
	//leo el fichero *_set.dat
	ifstream read(fichero_set);

	//si no lo encuentra o no puede habrirlo termina el programa
	if (!read){
		cout<<"No puedo abrir "<<fichero_set<<endl;
		exit(1);
	}

	//lee el contenido del archivo
	string cad="", sub_cad="";
	vector<float> fila;
	int size=0;//guarda el tamaño de la matriz
	while(!read.eof()) {
		++size;
		read >> cad;//leo por filas
		//recorro la fila caracter por caracter
		for(unsigned int i=0; i<cad.size(); ++i){
			if(cad.at(i)!=','){//si no encuentra una coma guarda el numero
				sub_cad+=cad.at(i);
			}else{//si encuentra una coma es que ya tiene el numero entero
				fila.push_back(stof(sub_cad));//lo añado al vector auxiliar
				sub_cad="";//y reseteo la sub cadena
			}
		}
		//añado el ultimo numero que no tiene delimitador a la derecha
		fila.push_back(stof(sub_cad));
		sub_cad="";//reseteo la subcadena

		//añado el vector de distancias de la fila leida
		atributos.push_back(fila);
		fila.clear();//limpio el vector
	}
	read.close();
	//borro el ultimo vector, ya que es la del salto de linea de la ultima linea, y contiene los mismos valores que la anterior
	atributos.pop_back();

	--size;//le decremento 1 porque cuenta el salto de linea del vacio
	//cout << size << endl;
	matriz.set_size(size,size);//redimensiono la matriz
	size_mat = size;

	read.open(fichero_set_const);//abro el archivo que contiene la matriz

	//si no lo encuentra o no puede habrirlo termina el programa
	if (!read){
		cout<<"No puedo abrir "<<fichero_set_const<<endl;
		exit(1);
	}

	//lee el contenido del archivo *_set_const
	cad="";
	string num="";
	int col=0, row=0, col2=1;
	while(!read.eof()) {
		read >> cad;
		//cout << row << ", " << col << endl;
		if(row<size){//esto lo hago porque lee el salto de linea de la ultima linea
			//recorro la linea
			for(unsigned int i=0; i<cad.size(); i+=2){
				//guardo el numero
				num = cad.at(i);
				if(num=="-"){//si es negativo guardo el - y el numero
					++i;
					num += cad.at(i);
				}
				//solo guardo la diagonal
				if(col>=col2){
					//guardo el valor en la matriz
					matriz(row,col)=stod(num);

					if(col != row and matriz(row,col) != 0){
						if(matriz(row,col) >0)
							ML.push_back(pair<int,int>(row,col));
						if(matriz(row,col) < 0)
							CL.push_back(pair<int,int>(row,col));
					}
				}

				++col;
			}
			//incremento de la fila e inicialización de la columna
			++row;
			col=0;
			col2=row+1;
		}

	}
	read.close();

	/*cout << "Resultados: " << endl;
	cout << matriz << endl;*/
}

//reset Centroides
void PAR::resetCentroides(){
	//vuelvo a generar aleatoriamente los centroides inicialmente distintos de dimensión n
	for(int i=0; i<k; ++i){
		centroides[i].clear();
		for(unsigned int e=0; e <atributos[0].size(); ++e)
			centroides[i].push_back(Rand());
	}
}

//imprime los atributos de cada nodo
void PAR::printDistanciasEuclideas(){
	for(vector<vector<float>>::iterator it=atributos.begin(); it != atributos.end(); ++it){
		for(vector<float>::iterator it2=(*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << endl;
	}
}
//imprime los centroides de cada cluster
void PAR::printCentroides(){
	int i=0;
	for(vector<vector<float>>::iterator it=centroides.begin(); it != centroides.end(); ++it){
		cout << i << " [ ";
		for(vector<float>::iterator it2=(*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << " ]" << endl;
		++i;
	}
}

//print the indexes of the attributes
void PAR::printRSI(){
	for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
		if(it+1 != RSI.end())
			cout << (*it) << ", ";
		else
			cout << (*it) << endl;
	}
}

//algoritmo Greedy
vector<vector<int>> PAR::algoritmoGreedy(){

	vector<vector<int>> cop_clusters(k);//vector de clusters copia para comparar con el cluster modificado
	vector<int> clusters_not_null;//vector of non-empty cluster indices
	vector<int> S(RSI.size(),-1);

	int pos=-1, iterations = 0, not_null=0;
	bool end = false;

	do{//mientras los vectores sean distintos
		if(clusters.size()>0){
			cop_clusters = clusters;//copio el vector de clusters antes de que sea modificado

			//clear vector of clusters
			clearClusters(false);

			clusters_not_null.clear();//clear the vector of clusters not empty
		}

		//go through all nodes
		for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){

			//save the cluster that best fits that node
			pos = minRestrictionsDistance(*it, clusters_not_null, S);
			//access to the cluster and add actual node
			clusters[pos].push_back(*it);
			S[*it] = pos;

			//if the actual cluster is not in the vector of not-null clusters
			if(find(clusters_not_null.begin(),clusters_not_null.end(),pos) == clusters_not_null.end())
				//add the cluster
				clusters_not_null.push_back(pos);

			pos = -1;//reset cluster
		}

		//update centroids
		for(int i = 0; i<k; ++i){
			//cout << "cluster: " << i << " size: " << clusters[i].size() << endl;
			if(clusters[i].size()>0){
				centroides[i] = updateDistance(clusters[i]);
			}
		}

		//cout << "calculate... " << iterations << endl;
		/*for(unsigned int i = 0; i< clusters.size(); ++i){
			cout << "cluster: " << i << endl;
			for(unsigned int e = 0; e < clusters[i].size(); ++e){
				if(clusters[i].size()>0)
					cout << clusters[i][e] << ", ";
			}
			cout << endl;
		}*/

		++iterations;
		//if the clusters don't undergo changes
		if(clusters==cop_clusters){

			//check that no cluster is empty
			for(int i =0; i < k; ++i){
				if(clusters[i].size()>0) ++not_null;
			}
			//if all the clusters aren't empty, finish Greedy
			if(not_null == k)
				end = true;
			//otherwise reset the centroids and clean clusters and cop_clusters
			else{
				resetCentroides();//reset the centroids
				not_null = 0;
				//clean clusters, clusters_not_null and cop_clusters
				clearClusters(false);
				for(unsigned int i=0; i< cop_clusters.size(); ++i){cop_clusters[i].clear();}
				clusters_not_null.clear();
			}
		}

	}while(not end);

	return clusters;//devuelvo el vector de cluster definitivo
}

//Update the distance
vector<float> PAR::updateDistance(vector<int> nodes){
	//save the actual distance
	vector<float> distance(atributos[nodes[0]].size(),0);
	for(unsigned int e=0; e<nodes.size(); ++e){
		//calculate average distance
		for(unsigned int i=0; i < atributos[nodes[e]].size(); ++i){
			//sumatorry
			distance[i] += atributos[nodes[e]][i];
		}
	}
	//average
	for(unsigned int i = 0; i < distance.size(); ++i){
		distance[i] = distance[i]/nodes.size();
	}

	return distance;
}

//calculate the closest and least restriction to cluster
int PAR::minRestrictionsDistance(int actual, vector<int> clusters_not_null, vector<int> S){
	int cluster=-1;
	float min_distance=999;//save the minimum distance and less restriction
	int less_restriction=999;
	float actual_distance=0;//save the actual distance
	int actual_restriction=0;//save the actual number of restrictions

	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){

		//calculate the Euclidea distance with the current cluster
		actual_distance = distanciaEuclidea(atributos[actual],centroides[i]);


		if(clusters_not_null.size()>0)
			actual_restriction = infeasibility(i, actual, clusters_not_null, S);

		//cout << "if( " << actual_restriction << " <= " << less_restriction << " and " << actual_distance << " <= " << min_distance << endl;
		//if the current cluster is less than the minimum saved, update the cluster and distance
		if(actual_restriction <=less_restriction){

			if((actual_restriction <less_restriction) || (actual_distance < min_distance && actual_restriction <=less_restriction)){

				min_distance = actual_distance;
				less_restriction = actual_restriction;
				cluster = i;
				//cout << "YES " << min_distance << ", " << less_restriction << endl;
			}
		}
	}
	/*if(less_restriction >0)
		cout << "MIN RESTRICTIONS: " << less_restriction << endl;*/
	return cluster;
}

//calculate the cluster with minimum distance, only first iteration
int PAR::minDistance(int actual){
	int cluster=-1;//save the cluster with minimum distance
	float min_distance=999;//save the minimum distance
	float actual_distance=0;//save the actual distance

	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){
		//calculate the Euclidea distance with the current cluster
		actual_distance = distanciaEuclidea(atributos[actual],centroides[i]);
		//cout << " distancia euclidea " << actual_distance << endl;

		//if the current cluster is less than the minimum saved, update the cluster and distance
		if(actual_distance < min_distance){
			min_distance = actual_distance;
			cluster = i;
		}
	}
	//cout << " minima distancia euclidea " << min_distance << " con cluster: " << cluster << endl;

	return cluster;//return the cluster with minimum distance
}

//calculate infeasibility when assigning an atribute to each cluster and return the minimum
int PAR::infeasibility(int clust, int actual, vector<int> clusters_not_null, vector<int> S){
	// number of restrictions, matrix column and row and not empty cluster indexes
	int rest=0, col=-1, row=-1;
	unsigned int max = CL.size();;
	if(ML.size()>CL.size())
		max = ML.size();

	for (unsigned int i = 0; i<max; ++i) {
		if(i < ML.size()){
			if(ML[i].first == actual){
				if(S[ML[i].second] != -1 && S[ML[i].second] != clust){
					++rest;
				}
			}else if(ML[i].second == actual){
				if(S[ML[i].first] != -1 && S[ML[i].first] != clust){
					++rest;
				}
			}
		}

		if(i < CL.size()){
			if(CL[i].first == actual){
				if(S[CL[i].second] != -1 && S[CL[i].second] == clust){
					++rest;
				}
			}else if(CL[i].second == actual){
				if(S[CL[i].first] != -1 && S[CL[i].first] == clust){
					++rest;
				}
			}
		}
	}

	//walk through non-empty clusters
	/*for(unsigned int i=0; i<clusters_not_null.size(); ++i){
		//walk through nodes of cluster
		for(unsigned int e = 0; e<clusters[clusters_not_null[i]].size(); ++e){

			//Is a triangular matrix
			if(clusters[clusters_not_null[i]][e]<actual){
				col = clusters[clusters_not_null[i]][e];
				row = actual;
			}else{
				col = actual;
				row = clusters[clusters_not_null[i]][e];
			}
			////if it isn't the identity and has a constraint
			if(matriz(col,row)!=0 and col != row){
				//if the node is not in the current cluster and the constraint is ML
				if(matriz(col,row) == 1 and (int) clusters_not_null[i] != clust)
					++rest;//increase constraint
				//if the node is in the current cluster and the constraint is CL
				else if(matriz(col,row) == -1 and (int) clusters_not_null[i] == clust)
					++rest;
			}
		}
	}*/

	return rest;// return the number of restrictions it violates
}

int PAR::infeasibility(vector<int> S){
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

//calcula la distancia euclidea entre 2 nodos
float PAR::distanciaEuclidea(vector<float> nod1, vector<float> nod2){
	//la formula es: sqrt(sumatoria((a_i - b_i)²))
	float suma=0;

	//sumatoria((a_i - b_i)²)
	for(unsigned int i=0; i<nod1.size();++i){
		suma += pow((nod1[i]-nod2[i]),2);
	}
	//sqrt(sumatoria)
	return sqrt(suma);
}

//random assignment of each node with a cluster
void PAR::randomAssign(){
	bool not_null = true;
	stack<int> clusters_null;

	//go through all nodes
	for(unsigned int i = 0; i < RSI.size(); ++i){
		//if it has traversed half the nodes
		if(i == RSI.size()/2){
			//check that no cluster is empty
			for(unsigned int i = 0; i < clusters.size(); ++i){
				//if the cluster is empty
				if(clusters[i].size() == 0){//add the cluster to the stack
					not_null = false;
					clusters_null.push(i);
				}
			}
		}


		//if all clusters aren't empty or haven't yet traveled half of the nodes
		if(not_null)
			clusters[rand() % k + 0].push_back(RSI[i]);//randomly assign a cluster
		else{//else it asigns to the cluster that is empty
			clusters[clusters_null.top()].push_back(RSI[i]);
			clusters_null.pop();
		}
		//if all the clusters aren't empty, assign the rest of the nodes randomly
		if(clusters_null.size() == 0)
			not_null = true;

	}
}

void PAR::clearClusters(bool all){

	for(unsigned int i = 0; i<clusters.size(); ++i){
		clusters[i].clear();
	}
	if(all)
		clusters.clear();
}

//create the vector of clusters assigned to each node
vector<int> PAR::createS(){
	vector<int> S(RSI.size());

	for(unsigned int e = 0; e < clusters.size(); ++e){
		for(vector<int>::iterator it = clusters[e].begin(); it != clusters[e].end(); ++it)
			S[*it] = e;
	}
	return S;
}

vector<vector<int>> PAR::algoritmoBL(){

	vector<int> S = createS();//create the vector S
	int infease = infeasibility(S);//calculate infeasibility
	int gen_deviation = generalDeviation(clusters); //calculate General Deviation
	float landa = createLanda();//calculate landa
	//and calculate the fitness
	float f = gen_deviation + (infease * landa);
	float new_f;

	//copy of S and clusters
	vector<int> S_cop = S;
	vector<vector<int>> clusters_cop = clusters;
	int i = 0;

	bool end = false;
	int iterate = 0;

	do{
		//if the node to change the cluster, the cluster has more than 1 element
		if(clusters_cop[S[i]].size()>0){
			//change of cluster
			S_cop[i] = rand() % k + 0;
			clusters_cop[S[i]].erase(find(clusters_cop[S[i]].begin(),clusters_cop[S[i]].end(),i));
			clusters_cop[S_cop[i]].push_back(i);

			//calculate the new fitness
			infease = infeasibility(S_cop);
			gen_deviation = generalDeviation(clusters_cop);

			new_f = gen_deviation + (infease * landa);

			//if the new f is less than the current f
			if(new_f < f){
				//update f and the clusters
				f = new_f;
				clusters = clusters_cop;
				S = S_cop;
			//if the new f is the same as the current f and there aren't change in the clusters
			}else if(new_f == f && clusters == clusters_cop && S == S_cop){
				end = true;//it end
			}else{
				clusters_cop = clusters;
				S_cop = S;
			}
		}
		++i;

	}while(!end || iterate < 100000);
	return clusters;
}

//calculate the general deviation
float PAR::generalDeviation(vector<vector<int>> v_clust){
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

//calculate Landa
float PAR::createLanda(){
	float lan = 0, actual_distance=0;
	//float node1=-1, node2=-1;
	//calculate the maximum distance
	for(unsigned int i = 0; i < atributos.size(); ++i){
		for(unsigned int e = i+1; e < atributos.size(); ++e){
			actual_distance = distanciaEuclidea(atributos[i],atributos[e]);
			if(actual_distance > lan){
				lan = actual_distance;
				/*node1 = i;
				node2 = e;*/
			}
		}
	}
	//count the number of restrictions total
	uvec rest = (find(matriz == 1 or matriz == -1));
	//cout << rest.size() << endl;

	//count number of restriction betwen node1 and node 2
	//cout << endl << node1 << " and " << node2 << ", distance: " << lan << endl;
	lan = lan/rest.size();

	//cout << "landa: " << lan << endl;

	return lan;
}
