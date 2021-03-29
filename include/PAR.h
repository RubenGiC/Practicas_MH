/*
 * lector_ficheros.h
 *
 *  Created on: 14 mar. 2021
 *      Author: Ruben Girela Castellón
 */
#ifndef INCLUDE_PAR_H_
#define INCLUDE_PAR_H_

#include <iostream>
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <armadillo>
using namespace std;
using namespace arma;

class PAR{

public:

	vector<vector<float>> atributos;//creo el vector de atributos
	mat matriz;//creo la matriz vacia de restricciones
	int k=7; //es el numero de cluster (es un numero fijo)
	vector<vector<float>> centroides;//vector de centroides
	vector<vector<int>> clusters;//es la asignación de cada elementoa cada cluster
	vector<int> RSI; //attribute index vector
	int size_mat;
	vector<int> clusters_not_null;//vector of non-empty cluster indices


	PAR(string fichero_set, string fichero_set_const);//constructor que inicializa los valores
	void lectura(string fichero_set, string fichero_set_const);//lee los archivos
	void print_distancias_euclideas();//imprime el vector de distancias euclideas
	void print_centroides();//imprime el vector de centroides
	void print_RSI();
	vector<vector<int>> algoritmo_greedy();//algoritmo greedy
	//calculate the closest and least restriction cluster
	int min_restrictions(int actual);
	//calculate the cluster with minimum distance, only in the first iteration
	int min_distance(int actual);


	//CALCULO DEL INFEASIBILITY Y DEVUELVO EL NUMERO DE RESTRICCIONES QUE INCUMPLE DEL CLUSTER QUE RECIBE
	int infeasibility(int cluster, int actual);;


	//CALCULO QUE NODO TIENE LA DISTANCIA EUCLIDEA MINIMA
	float distancia_euclidea(vector<float> nod1, vector<float> nod2);//calcula la distancia de 2 puntos
	float suma_total_distancias(int nodo);//suma todas las distancias de todos los nodos con la de un nodo concreto
	int min_distancia();//calcula que nodo tiene menor distancia

};
#endif /* INCLUDE_PAR_H_ */
