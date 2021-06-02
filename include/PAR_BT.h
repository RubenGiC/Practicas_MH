/*
 *
 *  Created on: 20 mar. 2021
 *      Author: Ruben Girela Castellón
 */
#ifndef INCLUDE_PAR_BT_H_
#define INCLUDE_PAR_BT_H_

#include <iostream>
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <armadillo>
using namespace std;
using namespace arma;

class PARBT{

public:

	vector<vector<float>> atributos;//creo el vector de atributos
	mat matriz;//creo la matriz vacia de restricciones
	int k=7; //es el numero de cluster (es un numero fijo)
	vector<vector<float>> centroides;//vector de centroides
	vector<int> RSI; //attribute index vector
	int size_mat;//size of matriz
	vector<pair<int,int>> CL;//save the par of nodes that it have CL constraints
	vector<pair<int,int>> ML;//save the par of nodes that it have ML constraints
	vector<int> S;//list of clusters assigned to each node
	float landa;

	PARBT(string fichero_set, string fichero_set_const, int semilla);//constructor que inicializa los valores
	void lectura(string fichero_set, string fichero_set_const);//lee los archivos

	//PRINT TO CONSOLE
	void printDistanciasEuclideas();//imprime el vector de distancias euclideas
	void printCentroides();//imprime el vector de centroides
	void printRSI();

	//ALGORITHMS
	vector<int> algoritmoBL(vector<int> S_cop, int max_iter, int &iterations, float &f);//BL algorithm
	vector<int> BMB(int max_iter, int n_solutions);//algoritmo de Busqueda Multiarranque Básica
	vector<int> ILS(int max_iter, int n_iterations);//algoritmo de Busqueda Local Reiterada
	vector<int> ES(int max_iter, float mu, float fi, float tf);//algoritmo de Enfriamiento Simulado
	vector<int> ILS_ES(int max_iter, float mu, float fi, float tf);//algoritmo Hibrido entre el Enfriamiento Simulado y el de Busqueda Local Reiterada

	//uniform mutation operator
	vector<int> fixedSegmentMutation(const vector<int> &sol);

	//MINIMIZATION FUNCTIONS
	//calculate the closest and least restriction cluster
	int minRestrictionsDistance(int actual, bool first);

	//MAXIMIZATION FUNCTION
	//calculate Landa
	float createLanda();
	//calculate general deviation
	float generalDeviation(const vector<int> &s_cop);

	//CALCULO DEL INFEASIBILITY Y DEVUELVO EL NUMERO DE RESTRICCIONES QUE INCUMPLE DEL CLUSTER QUE RECIBE
	int infeasibility(int cluster, int actual);
	int infeasibility(vector<int> S_cop);//return the number of restrictions that the solution has

	//CALCULO QUE NODO TIENE LA DISTANCIA EUCLIDEA MINIMA
	float distanciaEuclidea(vector<float> nod1, vector<float> nod2);//calcula la distancia de 2 puntos

	//generates the possible neighborhoods
	vector<pair<int,int>> generateNeig(vector<int> S_cop);

	//calculate a better fitness
	vector<int> betterFitness(const vector<int> &S_cop, const vector<pair<int,int>> &vecindario, float &f, float landa, int &it, int max);

	//calculate the fitness
	float fitness(const vector<int> &solution);

	//OTHER FUNCTIONS
	//reset centroides
	void resetCentroides();
	//random assignment solution
	vector<vector<int>> randomAssign(int n);
	//genearate random solution
	vector<int> randomSolution();
	//update the distance for each cluster
	vector<vector<float>> updateDistance(const vector<int> &nodes);
	//find all elements of the cluster
	vector<int> findInCluster(const vector<int> &s_cop, int clust);
	//count all elements of the cluster
	int countCluster(const vector<int> &s_cop, int clust);

	void shuffleRSI();
};
#endif /* INCLUDE_PAR_BT_H_ */
