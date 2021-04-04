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

	//PRINT TO CONSOLE
	void printDistanciasEuclideas();//imprime el vector de distancias euclideas
	void printCentroides();//imprime el vector de centroides
	void printRSI();

	//ALGORITHMS
	vector<vector<int>> algoritmoBL();//BL algorithm
	vector<vector<int>> algoritmoGreedy();//algoritmo greedy

	//MINIMIZATION FUNCTIONS
	//calculate the closest and least restriction cluster
	int minRestrictionsDistance(int actual);
	//calculate the cluster with minimum distance, only in the first iteration
	int minDistance(int actual);

	//MAXIMIZATION FUNCTION
	//calculate Landa
	float landa();

	//CALCULO DEL INFEASIBILITY Y DEVUELVO EL NUMERO DE RESTRICCIONES QUE INCUMPLE DEL CLUSTER QUE RECIBE
	int infeasibility(int cluster, int actual);;

	//CALCULO QUE NODO TIENE LA DISTANCIA EUCLIDEA MINIMA
	float distanciaEuclidea(vector<float> nod1, vector<float> nod2);//calcula la distancia de 2 puntos

	//OTHER FUNCTIONS
	//update the distance for each cluster
	vector<float> updateDistance(vector<int> nodes);
	//reset centroides
	void resetCentroides();
	//random assignment
	void randomAssign();
	//clear the clusters
	void clearClusters(bool all=true);
};
#endif /* INCLUDE_PAR_H_ */
