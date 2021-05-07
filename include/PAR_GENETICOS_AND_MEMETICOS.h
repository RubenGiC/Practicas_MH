/*
 * lector_ficheros.h
 *
 *  Created on: 23 april 2021
 *      Author: Ruben Girela Castellón
 */
#ifndef INCLUDE_PAR_GENETICOS_AND_MEMETICOS_H_
#define INCLUDE_PAR_GENETICOS_AND_MEMETICOS_H_

#include <iostream>
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <armadillo>
using namespace std;
using namespace arma;

enum TIPE_CROSS{
	AGG_UN, AGG_SF, AGE_UN, AGE_SF
};

class PAR_GM{

public:

	vector<vector<float>> atributos;//creo el vector de atributos
	mat matriz;//creo la matriz vacia de restricciones
	int k=7; //es el numero de cluster (es un numero fijo)

	vector<int> RSI; //attribute index vector
	int size_mat;//size of matriz
	vector<pair<int,int>> CL;//save the par of nodes that it have CL constraints
	vector<pair<int,int>> ML;//save the par of nodes that it have ML constraints
	vector<int> S;//list of clusters assigned to each node
	int size_sol;
	float landa;

	vector<vector<int>> vector_solutions; //Vector of solutions

	PAR_GM(string fichero_set, string fichero_set_const, int semilla);//constructor que inicializa los valores
	void lectura(string fichero_set, string fichero_set_const);//lee los archivos

	//PRINT TO CONSOLE
	void printDistanciasEuclideas();//imprime el vector de distancias euclideas
	void printCentroides(vector<vector<float>> centroides);//imprime el vector de centroides
	void printRSI();
	void printS(int s=0);
	void printSolution(vector<int> s);

	//ALGORITHMS
	vector<vector <int>> AGG(TIPE_CROSS cruce, float probability, int stop);
	vector<vector <int>> AGE(TIPE_CROSS cruce, float probability, int stop);
	vector<int> GENETIC(TIPE_CROSS tipo, float probability, int stop);

	//use the binary tournament, to select the best
	vector<vector<int>> selectionOperator(vector<vector<int>> actual, int tourney);

	//CROSSOVER OPERATORS
	//uniform crossover operator
	vector<vector<int>> uniformCross(vector<vector<int>> padres, float probability);
	//fixed segment crossover operator
	vector<vector<int>> fixedSegmentCross(vector<vector<int>> padres, float probability);
	//uniform mutation operator
	vector<vector<int>> uniformMutation(vector<vector<int>> padres);
	//replacement operator
	vector<vector<int>> replaceOperator(vector<vector<int>> hijos);

	//MINIMIZATION FUNCTIONS
	//calculate the closest and least restriction cluster
	int minRestrictionsDistance(int actual, bool first);

	//MAXIMIZATION FUNCTION
	//calculate Landa
	float createLanda();
	//calculate general deviation
	float generalDeviation(vector<int> s_cop);

	//CALCULO DEL INFEASIBILITY Y DEVUELVO EL NUMERO DE RESTRICCIONES QUE INCUMPLE DEL CLUSTER QUE RECIBE
	int infeasibility(int cluster, int actual);
	int infeasibility(vector<int> S_cop);//return the number of restrictions that the solution has

	//CALCULO QUE NODO TIENE LA DISTANCIA EUCLIDEA MINIMA
	float distanciaEuclidea(vector<float> nod1, vector<float> nod2);//calcula la distancia de 2 puntos

	//check wich solution is better
	int betterFitness(vector<vector<int>> padres, int indv1, int indv2);
	//calculate fitness
	float fitness(vector<int> solution);

	//OTHER FUNCTIONS
	//update the distance for each cluster
	vector<vector<float>> updateDistance(vector<int> nodes);
	//find all elements of the cluster
	vector<int> findInCluster(vector<int> s_cop, int clust);
	//reset centroides
	//void resetCentroides();
	//random assignment
	void randomAssign(int n);
	//clear the clusters
	void clearClusters(bool all=true);
	void shuffleRSI();
};
#endif /* INCLUDE_PAR_H_ */
