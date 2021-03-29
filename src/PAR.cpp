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

	//insert all indces of atributos
	for(unsigned int i=0; i<atributos.size(); ++i){
		RSI.push_back(i);
	}

	//despues barajo los indices de los atributos
	srand(unsigned (37));//genero una semilla fija
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
				if(col>=col2)
					//guardo el valor en la matriz
					matriz(row,col)=stod(num);

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

//imprime los atributos de cada nodo
void PAR::print_distancias_euclideas(){
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
void PAR::print_centroides(){
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
void PAR::print_RSI(){
	for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
		if(it+1 != RSI.end())
			cout << (*it) << ", ";
		else
			cout << (*it) << endl;
	}
}

//algoritmo Greedy
vector<vector<int>> PAR::algoritmo_greedy(){

	vector<vector<int>> cop_clusters;//vector de clusters copia para comparar con el cluster modificado
	int pos=-1;

	do{//mientras los vectores sean distintos
		if(clusters.size()>0)//si el vector de clusters no esta vacio
			cop_clusters = clusters;//copio el vector de clusters antes de que sea modificado
		for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
			if(it != RSI.begin()){
				pos = min_restrictions(*it);
				clusters[pos].push_back(*it);
				if(find(clusters_not_null.begin(),clusters_not_null.end(),pos) != clusters_not_null.end())
					clusters_not_null.push_back(pos);
			}else{
				pos = min_distance(*it);
				clusters[pos].push_back(*it);
				clusters_not_null.push_back(pos);
			}
			pos = -1;
		}

	}while(clusters!=cop_clusters);

	return clusters;//devuelvo el vector de cluster definitivo
}

//calculate the closest and least restriction to cluster
int PAR::min_restrictions(int actual){
	int cluster=-1;
	float min_distance=999;//save the minimum distance and less restriction
	int less_restriction=999;
	float actual_distance=0;//save the actual distance
	int actual_restriction=0;//save the actual number of restrictions
	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){
		//calculate the Euclidea distance with the current cluster
		actual_distance = distancia_euclidea(atributos[actual],centroides[i]);

		actual_restriction = infeasibility(i, actual);

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
	//cout << min_distance << " <-> " << less_restriction << " <-> " << endl;
	return cluster;
}

//calculate the cluster with minimum distance, only first iteration
int PAR::min_distance(int actual){
	int cluster=-1;//save the cluster with minimum distance
	float min_distance=999;//save the minimum distance
	float actual_distance=0;//save the actual distance

	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){
		//calculate the Euclidea distance with the current cluster
		actual_distance = distancia_euclidea(atributos[actual],centroides[i]);
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
int PAR::infeasibility(int clust, int actual){
	// number of restrictions, matrix column and row and not empty cluster indexes
	int rest=0, col=-1, row=-1, pos=-1;

	//vector that save the node and constraint
	vector<pair<int,int>> restriction_clust;

	//-1 CL (Cannot-Link) y 1 ML (Must-Link)

	//traverses the row of the constraint matrix of the current node
	for(int i =0; i<size_mat; ++i){

		//Is a triangular matrix
		if(i<actual){
			col = i;
			row = actual;
		}else{
			col = actual;
			row = i;
		}

		//if it isn't the diagonal and has a constraint
		if(col != row && matriz(i,actual)!=0){
			//walk through non-empty clusters
			for(unsigned int e=0; e<clusters_not_null.size(); ++e){
				//saves the current position of the non-emplty cluster
				pos = clusters_not_null[e];

				//if found by the node in current cluster
				if(find(clusters[pos].begin(),clusters[pos].end(),i)!=clusters[pos].end()){
					//saves node and constraint
					restriction_clust.push_back(make_pair(i,matriz(col,row)));
				}
			}
		}
	}
	//if the cluster is non-empty
	if(clusters[clust].size()>0){

		//walk through the vector with the constraints
		for(unsigned int i=0; i< restriction_clust.size(); ++i){

			//if the node is in the current cluster
			if(find(clusters[clust].begin(),clusters[clust].end(),restriction_clust[i].first) != clusters[clust].end()){
				//and the constraint is CL
				if(restriction_clust[i].second == -1)
					++rest;//increased constraint
			}else{//in case it isn't and the constraint is ML
				if(restriction_clust[i].second == 1)
					++rest;//increased constraint
			}
		}
	}else{//if the cluster is empty
		//walk through the vector with the constraints
		for(unsigned int i=0; i< restriction_clust.size(); ++i){
			//and increase the number of constraints if the constraint is ML
			if(restriction_clust[i].second == 1)
				++rest;
		}
	}
	return rest;// return the number of restriction
}

//calcula la distancia euclidea entre 2 nodos
float PAR::distancia_euclidea(vector<float> nod1, vector<float> nod2){
	//la formula es: sqrt(sumatoria((a_i - b_i)²))
	float suma=0;

	//sumatoria((a_i - b_i)²)
	for(unsigned int i=0; i<nod1.size();++i){
		suma += pow((nod1[i]-nod2[i]),2);
	}
	//sqrt(sumatoria)
	return sqrt(suma);
}

//suma todas las distancias euclideas de un nodo concreto
float PAR::suma_total_distancias(int nodo){
	float suma_total=0;
	//recorro todos los nodos menos si mismo
	for(unsigned int i=0;i<atributos.size(); ++i){
		if(i!=(unsigned int) nodo){
			//calculo la distancia euclidea entre ellos y la sumo
			suma_total += distancia_euclidea(atributos[nodo],atributos[i]);
		}
	}
	return suma_total;
}
//obtiene el nodo con menor distancia
int PAR::min_distancia(){
	int pos_nod=-1;//guardo el nodo con menor distancia
	float dist=0, dist_min=999999;//guardo la distancia actual y la distancia minima

	for(unsigned int i=0;i<atributos.size(); ++i){

		//calcula la distancia euclidea de todos los nodos con respecto al nodo actual
		dist = suma_total_distancias(i);
		//cout << "nodo: " << i << ", distancia: " << dist << endl;
		//si la distancia actual es mas pequeña que la minima encontrada
		if(dist < dist_min){
			//actualizo la distancia y el nodo
			dist_min = dist;
			pos_nod = i;
		}
	}
	cout << "El nodo: " << pos_nod << " tiene la distancia minima de " << dist_min << endl;
	return pos_nod;
}
