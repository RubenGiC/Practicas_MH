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
	/*PAR *par_zoo20 = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_20.const");
	PAR *par_glass10 = new PAR("datos/glass_set.dat", "datos/glass_set_const_10.const");
	PAR *par_glass20 = new PAR("datos/glass_set.dat", "datos/glass_set_const_20.const");
	PAR *par_bupa10 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_10.const");
	PAR *par_bupa20 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_20.const");*/
	vector<vector<int>> clusters_sol;
	ResultAlgorithms *results = new ResultAlgorithms();

	string results_all;

	for(int i = 0; i<5; ++i){

		clock_t start = clock();
		clusters_sol = par_zoo10->algoritmoGreedy();
		clock_t end = clock();
		float elapsed = float(end - start)/CLOCKS_PER_SEC;
		results_all += to_string(i+1) + " Elapsed (Greedy PAR zoo 10): " + to_string(elapsed) + "(seconds)\n";
		cout << (i+1) << " Elapsed (Greedy PAR zoo 10): " << elapsed << "(seconds)" << endl;

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


		results_all += "\tInfeas: " + to_string(results->Infeasable(clusters_sol, par_zoo10->ML, par_zoo10->CL)) + "\n";
		results_all += "\t Error Distance: " + to_string(abs(results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides)-0.9048)) + "\n";
		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;
		par_zoo10->clearClusters(false);//clear the clusters
		par_zoo10->shuffleRSI();
		par_zoo10->resetCentroides();
	}
	cout << results_all << endl;

	//PRUEBAS


	//leo los archivos
	//"datos/bupa_set.dat", "datos/bupa_set_const_10.const"
	PAR *par = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_10.const");
	//par->lectura("datos/bupa_set.dat", "datos/bupa_set_const_10.const");

	//clusters_sol = par->algoritmoGreedy();

	//BL (LOCAL SEARCH)
	par->clearClusters(false);//clear the clusters
	par->resetCentroides();//randomly generate the centroids
	par->randomAssign();//assign each node to a cluster randomly

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
