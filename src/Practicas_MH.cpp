//============================================================================
// Name        : Practicas_MH.cpp
// Author      : Ruben Girela Castellón
// Version     :
// Copyright   : Ruben Girela Castellón
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <armadillo>
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <string.h>

#include "../include/PAR.h"
#include "../include/result_algorithms.h"
using namespace std;
using namespace arma;

int main() {
	PAR *par_zoo10 = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_10.const");
	//PAR *par_zoo20 = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_20.const");
	PAR *par_glass10 = new PAR("datos/glass_set.dat", "datos/glass_set_const_10.const");
	//PAR *par_glass20 = new PAR("datos/glass_set.dat", "datos/glass_set_const_20.const");
	PAR *par_bupa10 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_10.const");
	//PAR *par_bupa20 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_20.const");
	vector<vector<int>> clusters_sol;
	ResultAlgorithms *results = new ResultAlgorithms();

	string zoo_20 = "", zoo_10 = "", glass_10 = "", glass_20 = "", bupa_10 = "", bupa_20 = "";
	float elapsed;
	clock_t start;
	clock_t end;

	cout << "Calculate Greedy (ZOO, GLASS, BUPA) and BL (ZOO, GLASS, BUPA)" << endl;

	/*for(int i = 0; i<5; ++i){
		//ZOO 10
		start = clock();
		clusters_sol = par_zoo10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		zoo_10 += to_string(i+1) + " Elapsed (Greedy PAR ZOO 10): " + to_string(elapsed) + "(seconds)\n";

		zoo_10 += "\tInfeas: " + to_string(results->Infeasable(clusters_sol, par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";
		zoo_10 += "\tError Distance: " + to_string(abs(results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides)-0.904799856)) + "\n";
		zoo_10 += "\tFitness: " + to_string(results->Fitness(par_zoo10->atributos,par_zoo10->matriz, clusters_sol, par_zoo10->centroides, par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";
		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;

		cout << i << " ZOO ------------------------------" << endl;

		par_zoo10->clearClusters(false);//clear the clusters
		par_zoo10->shuffleRSI();
		par_zoo10->resetCentroides();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//GLASS 10
		start = clock();
		clusters_sol = par_glass10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		glass_10 += to_string(i+1) + " Elapsed (Greedy PAR GLASS 10): " + to_string(elapsed) + "(seconds)\n";

		glass_10 += "\tInfeas: " + to_string(results->Infeasable(clusters_sol, par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";
		glass_10 += "\tError Distance: " + to_string(abs(results->Distance(clusters_sol, par_glass10->atributos, par_glass10->centroides)-0.364290282)) + "\n";
		glass_10 += "\tFitness: " + to_string(results->Fitness(par_glass10->atributos,par_glass10->matriz, clusters_sol, par_glass10->centroides, par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";

		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;

		cout << i << " GLASS ------------------------------" << endl;

		par_glass10->clearClusters(false);//clear the clusters
		par_glass10->shuffleRSI();
		par_glass10->resetCentroides();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//BUPA 10
		start = clock();
		clusters_sol = par_bupa10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bupa_10 += to_string(i+1) + " Elapsed (Greedy PAR BUPA 10): " + to_string(elapsed) + "(seconds)\n";

		bupa_10 += "\tInfeas: " + to_string(results->Infeasable(clusters_sol, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";
		bupa_10 += "\tError Distance: " + to_string(abs(results->Distance(clusters_sol, par_bupa10->atributos, par_bupa10->centroides)-0.220423749)) + "\n";
		bupa_10 += "\tFitness: " + to_string(results->Fitness(par_bupa10->atributos,par_bupa10->matriz, clusters_sol, par_bupa10->centroides, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";

		par_bupa10->clearClusters(false);//clear the clusters
		par_bupa10->shuffleRSI();
		par_bupa10->resetCentroides();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();
		*/
		//ZOO 20
		/*start = clock();
		clusters_sol = par_zoo20->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		zoo_20 += to_string(i+1) + " Elapsed (Greedy PAR zoo 20): " + to_string(elapsed) + "(seconds)\n";

		zoo_20 += "\tInfeas: " + to_string(results->Infeasable(clusters_sol, par_zoo20->ML, par_zoo20->CL)) + "\n";
		zoo_20 += "\t Error Distance: " + to_string(abs(results->Distance(clusters_sol, par_zoo20->atributos, par_zoo20->centroides)-0.9048)) + "\n";
		zoo_20 += "\t Fitness: " + to_string(results->Fitness(par_zoo20->atributos,par_zoo20->matriz, clusters_sol, par_zoo20->centroides)) + "\n";
		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;
		par_zoo20->clearClusters(false);//clear the clusters
		par_zoo20->shuffleRSI();
		par_zoo20->resetCentroides();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();*/


		//cout << (i+1) << " Elapsed (Greedy PAR zoo 10): " << elapsed << "(seconds)" << endl;*/

		/*cout << "Solution clusters (Greedy PAR zoo 10):" << endl;
		int n = 0;
		for(vector<vector<int>>::iterator it = clusters_sol.begin(); it != clusters_sol.end(); ++it){
			cout << n << " [ ";
			for(vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
				if(it2+1 != (*it).end())
					cout << (*it2) << ", ";
				else
					cout << (*it2);
			}
			cout << " ]" << endl;
			++n;
		}*/
	//}

	cout << endl << "ZOO ********************************************* " << endl;
	cout << zoo_10 << endl;
	cout << endl << zoo_20 << endl;

	cout << "GLASS ********************************************* " << endl;

	cout << endl << glass_10 << endl;
	cout << endl << glass_20 << endl;

	cout << "BUPA ********************************************* " << endl;

	cout << endl << bupa_10 << endl;
	cout << endl << bupa_20 << endl;

	//PRUEBAS


	//leo los archivos
	//"datos/bupa_set.dat", "datos/bupa_set_const_10.const"
	PAR *par = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_10.const");
	//par->lectura("datos/bupa_set.dat", "datos/bupa_set_const_10.const");

	//BL (LOCAL SEARCH)
	par->clearClusters(false);//clear the clusters
	par->resetCentroides();//randomly generate the centroids
	par->randomAssign();//assign each node to a cluster randomly

	start = clock();
	clusters_sol = par->algoritmoBL();
	end = clock();
	elapsed = float(end - start)/CLOCKS_PER_SEC;
	cout << " Elapsed (BL PAR ZOO 10): " << elapsed << "(seconds)\n";

	cout << "\tInfeas: " << results->Infeasable(clusters_sol, par->ML, par->CL, par->createS()) << "\n";
	cout << "\tError Distance: " << abs(results->Distance(clusters_sol, par->atributos, par->centroides)-0.904799856) << "\n";
	cout << "\tFitness: " << results->Fitness(par->atributos,par->matriz, clusters_sol, par->centroides, par->ML, par->CL, par->createS()) << "\n";
	/*cout << "Asign cluster randomly:" << endl;
	n = 0;
	for(vector<vector<int>>::iterator it = par->clusters.begin(); it != par->clusters.end(); ++it){
		cout << n << " [ ";
		for(vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << " ]" << endl;
		++n;
	}*/
	return 0;
}
