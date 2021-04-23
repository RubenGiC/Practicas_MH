/*
 * lector_ficheros.cpp
 *
 *  Created on: 23 april 2021
 *      Author: Ruben Girela Castellón
 */
#include <iostream>
#include <armadillo>// para la matriz
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <math.h> // para calculos matematicos (sqrt, pow (x^y), ...)
#include <algorithm>// para barajar el vector
#include <ctime>// libreria de time para el shuffle aleatorio de vectores
#include <stack>// use stack

#include "../include/PAR_GENETICOS_AND_MEMETICOS.h"
#include "../include/random.h"
using namespace std;
using namespace arma;

PAR_GM::PAR_GM(string fichero_set, string fichero_set_const, int semilla){

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

	Set_random(semilla);//creo una semilla para los valores aleatorios

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
	srand(unsigned (semilla));//genero una semilla fija
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR_GM::shuffleRSI(){
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR_GM::lectura(string fichero_set, string fichero_set_const){
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
}

//reset Centroides
void PAR_GM::resetCentroides(){
	//vuelvo a generar aleatoriamente los centroides inicialmente distintos de dimensión n
	for(int i=0; i<k; ++i){
		centroides[i].clear();
		for(unsigned int e=0; e <atributos[0].size(); ++e)
			centroides[i].push_back(Rand());
	}
}

//print the elements of each cluster
void PAR_GM::printS(){
	vector<int> elements;
	int count = 0, total = 0;
	for(int i = 0; i<k; ++i){
		elements = findInCluster(S,i);
		cout << i << ": [ ";
		for(auto e = elements.begin(); e != elements.end(); ++e){
			cout << *e << ", ";
			++count;
		}
		cout << " ] n = " << count << endl;
		total += count;
		count = 0;
	}
	cout << "total elements = " << total << endl;
}

//imprime los atributos de cada nodo
void PAR_GM::printDistanciasEuclideas(){
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
void PAR_GM::printCentroides(){
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
void PAR_GM::printRSI(){
	for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
		if(it+1 != RSI.end())
			cout << (*it) << ", ";
		else
			cout << (*it) << endl;
	}
}

//algoritmo Greedy
vector<vector<int>> PAR_GM::algoritmoGreedy(){

	vector<vector<int>> cop_clusters(k);//vector de clusters copia para comparar con el cluster modificado
	vector<int> S_cop(RSI.size(),-1);
	S = S_cop;

	int pos=-1, not_null=0;
	//int iterations = 0;
	bool end = false, first = true;

	do{//mientras los vectores sean distintos
		if(clusters.size()>0){
			cop_clusters = clusters;//copio el vector de clusters antes de que sea modificado

			//clear vector of clusters
			clearClusters(false);
		}

		//go through all nodes
		for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
			//save the cluster that best fits that node
			pos = minRestrictionsDistance(*it, first);
			//access to the cluster and add actual node
			clusters[pos].push_back(*it);
			S[*it] = pos;

			pos = -1;//reset cluster
			first = false;
		}
		first = true;

		//update centroids
		for(int i = 0; i<k; ++i){
			//cout << "cluster: " << i << " size: " << clusters[i].size() << endl;
			if(clusters[i].size()>0){
				centroides[i] = updateDistance(clusters[i]);
			}
		}

		//++iterations;
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
				//clean clusters and cop_clusters
				clearClusters(false);
				for(unsigned int i=0; i< cop_clusters.size(); ++i){cop_clusters[i].clear();}
				//S = S_cop;
			}
		}

	}while(not end);

	return clusters;//devuelvo el vector de cluster definitivo
}

//Update the distance
vector<float> PAR_GM::updateDistance(vector<int> nodes){
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
int PAR_GM::minRestrictionsDistance(int actual, bool first){
	int cluster=-1;
	float min_distance=999;//save the minimum distance and less restriction
	int less_restriction=999;
	float actual_distance=0;//save the actual distance
	int actual_restriction=0;//save the actual number of restrictions

	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){

		//calculate the Euclidea distance with the current cluster
		actual_distance = distanciaEuclidea(atributos[actual],centroides[i]);

		//if it isn't the first node to enter
		if(!first)
			//calculates the number of constraints it violates
			actual_restriction = infeasibility(i, actual);

		//if the current cluster is less than the minimum saved, update the cluster and distance
		if(actual_restriction <=less_restriction){

			if((actual_restriction <less_restriction) || (actual_distance < min_distance && actual_restriction <=less_restriction)){

				min_distance = actual_distance;
				less_restriction = actual_restriction;
				cluster = i;
			}
		}
	}
	return cluster;
}

//calculate infeasibility when assigning an atribute to each cluster and return the minimum
int PAR_GM::infeasibility(int clust, int actual){
	// number of restrictions, matrix column and row and not empty cluster indexes
	int rest=0;
	//walk through each constrain ML
	for (unsigned int i = 0; i<ML.size(); ++i) {
		if(i < ML.size()){
			//if it find the current node in the constrain
			if(ML[i].first == actual){
				//and the second node has an assigned cluster and isn't in the same cluster
				if(S[ML[i].second] != -1 && S[ML[i].second] != clust){
					++rest;//increases the number of restrictions violated
				}
			}else if(ML[i].second == actual){
				if(S[ML[i].first] != -1 && S[ML[i].first] != clust){
					++rest;
				}
			}
		}
	}
	//walk through each constrain CL
	for (unsigned int i = 0; i<CL.size(); ++i) {
		//if it find the current node in the constrain
		if(CL[i].first == actual){
			//and the second node has an assigned cluster and is in the same cluster
			if(S[CL[i].second] != -1 && S[CL[i].second] == clust){
				++rest;//increases the number of restrictions violated
			}
		}else if(CL[i].second == actual){
			if(S[CL[i].first] != -1 && S[CL[i].first] == clust){
				++rest;
			}
		}
	}

	return rest;// return the number of restrictions it violates
}

//same the infeasibility(int clust, int actual) except it receives the solution
int PAR_GM::infeasibility(vector<int> S_cop){
	int restrictions = 0;
	//CALCULATE the max size
	unsigned int max = CL.size();
	if(ML.size()>max) max = ML.size();

	//walk through all restrictions
	for(unsigned int i = 0; i< max; ++i){
		if(i<CL.size()){
			//check if it violate the CL constraint
			if(S_cop[CL[i].first] == S_cop[CL[i].second])
				++restrictions;
		}

		if(i<ML.size()){
			//check if it violate the ML constraint
			if(S_cop[ML[i].first] != S_cop[ML[i].second])
				++restrictions;
		}
	}

	return restrictions;
}

//calcula la distancia euclidea entre 2 nodos
float PAR_GM::distanciaEuclidea(vector<float> nod1, vector<float> nod2){
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
void PAR_GM::randomAssign(){
	bool not_null = true;
	bool find = false;
	stack<int> clusters_null;

	S.resize(RSI.size());
	cout << S.size() << endl;

	//go through all nodes
	for(unsigned int i = 0; i < RSI.size(); ++i){
		//if it has traversed half the nodes
		if(i == RSI.size()/2){
			//check that no cluster is empty
			for(int e = 0; e < k; ++e){
				for(unsigned int j=0; j<RSI.size()/2 && !find; ++j){
					if(S[(int)RSI[j]] == e)
						find=true;
				}
				//if the cluster is empty
				if(!find){//add the cluster to the stack
					not_null = false;
					clusters_null.push(e);
				}else
					find = false;
			}
		}


		//if all clusters aren't empty or haven't yet traveled half of the nodes
		if(not_null)
			S[RSI[i]] = rand() % k + 0;
		else{//else it asigns to the cluster that is empty
			S[RSI[i]] = clusters_null.top();
			clusters_null.pop();
		}
		//if all the clusters aren't empty, assign the rest of the nodes randomly
		if(clusters_null.size() == 0)
			not_null = true;

	}
	//calculates centroids of that random solution
	for(int i = 0; i < k; ++i){
		centroides[i] = updateDistance(findInCluster(S,i));
	}
}

void PAR_GM::clearClusters(bool all){

	for(unsigned int i = 0; i<clusters.size(); ++i){
		clusters[i].clear();
	}
	if(all)
		clusters.clear();
}

//create the vector of clusters assigned to each node
vector<int> PAR_GM::createS(){
	S.resize(RSI.size());

	for(unsigned int e = 0; e < clusters.size(); ++e){
		for(vector<int>::iterator it = clusters[e].begin(); it != clusters[e].end(); ++it)
			S[*it] = e;
	}
	return S;
}

vector<vector<int>> PAR_GM::algoritmoBL(){
	int max = 100000;

	S = createS();//create the vector S

	int infease = infeasibility(S);//calculate infeasibility
	float gen_deviation = generalDeviation(S); //calculate General Deviation
	float landa = createLanda();//calculate landa
	//and calculate the fitness
	float f = gen_deviation + (infease * landa);

	vector<pair<int,int>> vecindario = generateNeig();//create the neighborhood
	random_shuffle(vecindario.begin(), vecindario.end());//shuffle vector

	bool end = false;//algorithm completion
	int iterate = 0;

	do{
		//if the number of iterate don't reach 100000
		if(iterate < max){//calculate the new solution
			iterate = betterFitness(vecindario, f, landa, iterate, max);
			//update centroides
			for(unsigned int i = 0; i < clusters.size(); ++i){
				centroides[i] = updateDistance(clusters[i]);
			}

			//calculate the new neighborhood
			vecindario = generateNeig();
			random_shuffle(vecindario.begin(), vecindario.end());//barajo el vector
		}else
			end = true;//else the algorithm ends

	}while(!end);
	//cout << iterate << endl;
	return clusters;
}

//genereate the possible neighborhoods
vector<pair<int,int>> PAR_GM::generateNeig(){
	vector<pair<int,int>> vecindario;
	for(unsigned int i=0; i<RSI.size(); ++i){
		for(int e = 0; e<k; ++e){
			vecindario.push_back(pair<int, int>(RSI[i],e));
		}
	}
	return vecindario;
}

//calculates the new better solution than the current solution
int PAR_GM::betterFitness(vector<pair<int,int>> vecindario, float &f, float landa, int it, int max){
	int iterate = it;
	int infease = 0;//calculate infeasibility
	float gen_deviation = 0; //calculate General Deviation
	//copy the solution and neighborhood
	vector<int> S_cop = S;
	vector<vector<int>> clusters_cop = clusters;
	float new_f = 0;

	//tour the neighborhood
	for(auto i=vecindario.begin(); i<vecindario.end(); ++i){
		//if the cluster has more than 1 node
		if(clusters_cop[S_cop[(*i).first]].size()>1){
			//and the new cluster is different from the current one
			if(S_cop[(*i).first] != (*i).second){
				//change
				S_cop[(*i).first] = (*i).second;

				//update clusters_cop
				//erase the node in old cluster
				clusters_cop[S[(*i).first]].erase(find(clusters_cop[S[(*i).first]].begin(),clusters_cop[S[(*i).first]].end(),(*i).first));
				//add the node in new cluster
				clusters_cop[S_cop[(*i).first]].push_back((*i).first);

				//calculate the new fitness
				gen_deviation = generalDeviation(S_cop);
				infease = infeasibility(S_cop);
				new_f = gen_deviation + (infease * landa);

				++iterate;//increase the number of iterations

				//if the new solution is better than current solution
				if(new_f < f){
					//update the current solution
					f = new_f;
					S = S_cop;
					clusters = clusters_cop;

					return iterate;//and return the actual iteration
				}
				//else restore to the previous solution
				S_cop = S;
				clusters_cop = clusters;

				//if the iterate reaches 100000
				if(iterate == max){
					return iterate;//ends
				}


			}
		}
	}
	//if it calculate all the possibilities, it ends
	return max;
}

//calculate the general deviation
float PAR_GM::generalDeviation(vector<int> s_cop){
	float distance = 0, intra_cluster = 0;
	vector<int> elements;

		//walk through each cluster
		for(int i = 0; i< k; ++i){

			elements = findInCluster(s_cop,i);
			//sumatorry(euclidean distance of all nodes in the cluster)
			for(vector<int>::iterator it = elements.begin(); it != elements.end(); ++it){
				distance += distanciaEuclidea(atributos[(*it)],centroides[i]);
			}
			//mean intra-cluster distance
			distance = distance / elements.size();
			intra_cluster += distance;
			distance = 0;
		}

		return intra_cluster/k;
}

//find all elements of the cluster
vector<int> PAR_GM::findInCluster(vector<int> s_cop, int clust){
	vector<int> elements;
	for(unsigned int i=0; i<s_cop.size(); ++i)
		if(clust == s_cop[i])
			elements.push_back(i);

	return elements;
}

//calculate Landa
float PAR_GM::createLanda(){
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
