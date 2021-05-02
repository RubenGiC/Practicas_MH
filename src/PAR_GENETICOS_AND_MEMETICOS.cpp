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
/*void PAR_GM::resetCentroides(){
	//vuelvo a generar aleatoriamente los centroides inicialmente distintos de dimensión n
	for(int i=0; i<k; ++i){
		centroides[i].clear();
		for(unsigned int e=0; e <atributos[0].size(); ++e)
			centroides[i].push_back(Rand());
	}
}*/

//print the elements of each cluster
void PAR_GM::printS(int s){
	vector<int> elements;
	int count = 0, total = 0;
	if(s == 0){

		for(int j=0; j< size_sol; ++j){
			cout << "SOLUCIÓN: " << j+1 << endl;
			for(int i = 0; i<k; ++i){
				elements = findInCluster(vector_solutions[j],i);
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
			total = 0;
		}
	}else{
		for(int i = 0; i<k; ++i){
			elements = findInCluster(vector_solutions[s],i);
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
void PAR_GM::printCentroides(vector<vector<float>> centroides){
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

vector<int> PAR_GM::GENETIC(TIPE_CROSS tipo){

	//choose the tipe of crossover operator and genetic algorithm
	if(tipo == (AGG_UN || AGG_SF))
		AGG(tipo);
	else
		AGE(tipo);



	return S;
}

vector <int> PAR_GM::AGE(TIPE_CROSS cruce){

	vector<vector<int>> vector_padres = selectionOperator(vector_solutions, 2);

	if(cruce == AGE_UN)
		uniformCross(vector_padres);
	else
		fixedSegmentCross(vector_padres);

	return S;
}

vector <int> PAR_GM::AGG(TIPE_CROSS cruce){

	vector<vector<int>> vector_padres = selectionOperator(vector_solutions, vector_solutions.size());

	if(cruce == AGG_UN)
		uniformCross(vector_padres);
	else
		fixedSegmentCross(vector_padres);

	return S;
}

vector<vector<int>> PAR_GM::selectionOperator(vector<vector<int>> actual, int tourney){
	int indv1= -1, indv2= -1;
	vector<vector<int>> padres;
	landa = createLanda();

	for(int i=0; i < tourney; ++i){
		do{
			indv1 = rand() % actual.size() + 0;
			indv2 = rand() % actual.size() + 0;
		}while(indv1 == indv2);

		//cout << indv1 << " vs " << indv2 << endl;

		cout << "THE BEST: " << betterFitness(actual,indv1,indv2) << endl;
		padres.push_back(actual[betterFitness(actual,indv1,indv2)]);
	}

	return padres;
}

//Update the distance
vector<vector<float>> PAR_GM::updateDistance(vector<int> nodes){
	//save the actual distance
	vector<float> distance(atributos[nodes[0]].size(),0);
	vector<vector<float>> centroides;
	centroides.resize(k);

	int size = 0;
	for(int j = 0; j<k; ++j){
		for(unsigned int e=0; e<nodes.size(); ++e){
			if(nodes[e] == j){
				++size;
				//calculate average distance
				for(unsigned int i=0; i < atributos[e].size(); ++i){
					//sumatorry
					distance[i] += atributos[e][i];
				}
			}
		}

		//average
		for(unsigned int i = 0; i < distance.size(); ++i){
			distance[i] = distance[i]/nodes.size();
		}
		centroides[j] = distance;
	}


	return centroides;
}

//calculate the closest and least restriction to cluster
int PAR_GM::minRestrictionsDistance(int actual, bool first){
	int cluster=-1;
	/*float min_distance=999;//save the minimum distance and less restriction
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
	}*/
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
void PAR_GM::randomAssign(int n){
	bool not_null = true;
	//bool find = false;
	stack<int> clusters_null;
	vector<int> clusters_asign;

	size_sol = n;
	vector_solutions.resize(n);

	S.resize(RSI.size());
	//cout << S.size() << endl;

	//initialize vector
	for(int i = 0; i<n; ++i){
		vector_solutions[i].resize(RSI.size());
	}

	for(int j=0; j<n; ++j){

		//go through all nodes
		for(unsigned int i = 0; i < RSI.size(); ++i){
			//if it has traversed half the nodes
			if(i == RSI.size()/2){
				//check that no cluster is empty
				for(int e = 0; e < k; ++e){
					/*for(unsigned int l=0; l<RSI.size()/2 && !find; ++l){
						if(vector_solutions[j][(int)RSI[l]] == e)
							find=true;
					}*/
					//if the cluster is empty
					/*if(!find){//add the cluster to the stack
						not_null = false;
						clusters_null.push(e);
					}else
						find = false;*/

					//find cluster empty
					if(find(clusters_asign.begin(), clusters_asign.end(), e) == clusters_asign.end()){

						clusters_null.push(e);
						not_null = false;
					}
				}
			}


			//if all clusters aren't empty or haven't yet traveled half of the nodes
			if(not_null){
				vector_solutions[j][RSI[i]] = rand() % k + 0;

				if(find(clusters_asign.begin(), clusters_asign.end(), vector_solutions[j][RSI[i]]) == clusters_asign.end()){
					clusters_asign.push_back(vector_solutions[j][RSI[i]]);
				}

			}else{//else it asigns to the cluster that is empty
				vector_solutions[j][RSI[i]] = clusters_null.top();
				clusters_null.pop();
			}
			//if all the clusters aren't empty, assign the rest of the nodes randomly
			if(clusters_null.size() == 0)
				not_null = true;

		}
	}
}

//calculates wich of the 2 individuals is the best
int PAR_GM::betterFitness(vector<vector<int>> padres, int indv1, int indv2){

	int infease = 0;//calculate infeasibility
	float gen_deviation = 0; //calculate General Deviation
	float f1 = 0, f2 = 0;

	//calculate the fitness indv1
	gen_deviation = generalDeviation(padres[indv1]);
	infease = infeasibility(padres[indv1]);
	f1 = gen_deviation + (infease * landa);

	//calculate the fitness indv2
	gen_deviation = generalDeviation(padres[indv2]);
	infease = infeasibility(padres[indv2]);
	f2 = gen_deviation + (infease * landa);

	cout << "Fitness: " << f1 << " vs " << f2 << endl;

	//if the new solution is better than current solution
	if(f1 < f2)
		return indv1;//and return the actual iteration

	return indv2;
}

//calculate the general deviation
float PAR_GM::generalDeviation(vector<int> s_cop){
	float distance = 0, intra_cluster = 0;
	vector<int> elements;
	vector<vector<float>> centroides;
	centroides = updateDistance(s_cop);

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
