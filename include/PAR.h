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
	int size_mat;//size of matriz
	vector<pair<int,int>> CL;//save the par of nodes that it have CL constraints
	vector<pair<int,int>> ML;//save the par of nodes that it have ML constraints
	vector<int> S;//list of clusters assigned to each node

	PAR(string fichero_set, string fichero_set_const, int semilla);//constructor que inicializa los valores
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
	int minRestrictionsDistance(int actual, bool first);

	//MAXIMIZATION FUNCTION
	//calculate Landa
	float createLanda();
	//calculate general deviation
	float generalDeviation(vector<vector<int>> v_clust);

	//CALCULO DEL INFEASIBILITY Y DEVUELVO EL NUMERO DE RESTRICCIONES QUE INCUMPLE DEL CLUSTER QUE RECIBE
	int infeasibility(int cluster, int actual);
	int infeasibility(vector<int> S_cop);//return the number of restrictions that the solution has

	//CALCULO QUE NODO TIENE LA DISTANCIA EUCLIDEA MINIMA
	float distanciaEuclidea(vector<float> nod1, vector<float> nod2);//calcula la distancia de 2 puntos

	//generates the possible neighborhoods
	vector<pair<int,int>> generateNeig();

	//calculate a better fitness
	int betterFitness(vector<pair<int,int>> vecindario, float &f, float landa, int it, int max);

	//OTHER FUNCTIONS
	//update the distance for each cluster
	vector<float> updateDistance(vector<int> nodes);
	//reset centroides
	void resetCentroides();
	//random assignment
	void randomAssign();
	//clear the clusters
	void clearClusters(bool all=true);
	void shuffleRSI();
	//create the vector of clusters assigned to each node
	vector<int> createS();
};
#endif /* INCLUDE_PAR_H_ */
