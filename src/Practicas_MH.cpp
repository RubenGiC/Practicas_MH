//============================================================================
// Name        : Practicas_MH.cpp
// Author      : Ruben Girela Castellón
// Version     :
// Copyright   : Ruben Girela Castellón
// Description : Apply Greedy algorithms and Local Search in the PAR problem (Grouping with Constraints Problem (GCP)) in C++, Ansi-style
//============================================================================

#include <iostream>
#include <armadillo>
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <string.h>

#include "../include/PAR.h"
#include "../include/PAR_GENETICOS_AND_MEMETICOS.h"
#include "../include/result_algorithms.h"

using namespace std;
using namespace arma;

int main(int argc, char** argv) {

	//int seed = 37;
	vector<int> seeds;
	vector<string> paths_data;
	vector<string> paths_const;
	int iterations = 100000;

	string tipo;

	if(argc>1){
		if(argc > 2){
			if(argc%2 == 0){
				for(int i = 1; i <argc-1; i+=3){
					//seed
					seeds.push_back(atoi(argv[i]));
					//data
					paths_data.push_back(argv[i+1]);
					//constraints
					paths_const.push_back(argv[i+2]);
				}
			}else{
				for(int i = 1; i <argc-1; i+=3){
					//seed
					seeds.push_back(atoi(argv[i]));
					//data
					paths_data.push_back(argv[i+1]);
					//constraints
					paths_const.push_back(argv[i+2]);
				}
			}
		}
		if(argc%2 == 0)
			iterations = atoi(argv[argc-1]);
	}

	/*for(unsigned int i = 0; i<paths_data.size(); ++i){
		cout << "SEED: " << seeds[i] << endl;
		cout << "PATH DATA: " << paths_data[i] << endl;
		cout << "PATH CONST: " << paths_const[i] << endl;
	}*/

	/*PAR *par_zoo10 = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_10.const", seed);
	PAR *par_zoo20 = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_20.const", seed);
	PAR *par_glass10 = new PAR("datos/glass_set.dat", "datos/glass_set_const_10.const", seed);
	PAR *par_glass20 = new PAR("datos/glass_set.dat", "datos/glass_set_const_20.const", seed);
	PAR *par_bupa10 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_10.const", seed);
	PAR *par_bupa20 = new PAR("datos/bupa_set.dat", "datos/bupa_set_const_20.const", seed);
	vector<vector<int>> clusters_sol, clusters_sol_gred, clusters_sol_bl;
	ResultAlgorithms *results = new ResultAlgorithms();

	string greedy_zoo_20 = "", greedy_zoo_10 = "", greedy_glass_10 = "", greedy_glass_20 = "", greedy_bupa_10 = "", greedy_bupa_20 = "";
	string bl_zoo_20 = "", bl_zoo_10 = "", bl_glass_10 = "", bl_glass_20 = "", bl_bupa_10 = "", bl_bupa_20 = "";
	float elapsed;
	clock_t start;
	clock_t end;
	clock_t start_global;
	clock_t end_global;

	cout << "Calculate Greedy (ZOO, GLASS, BUPA) Aproximate 1 minute" << endl;
	start_global = clock();
	for(int i = 0; i<5; ++i){
		//ZOO 10
		start = clock();
		clusters_sol = par_zoo10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_zoo_10 += to_string(i+1) + " Elapsed (Greedy PAR ZOO 10): " + to_string(elapsed) + "(seconds)\n";

		greedy_zoo_10 += "\tInfeas: " + to_string(results->Infeasable(par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";
		greedy_zoo_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides,"ZOO")) + "\n";
		greedy_zoo_10 += "\tFitness: " + to_string(results->Fitness(par_zoo10->atributos,par_zoo10->matriz, clusters_sol, par_zoo10->centroides, par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";
		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;

		//cout << i << " ZOO GREEDY 10 ------------------------------" << endl;

		par_zoo10->clearClusters(false);//clear the clusters
		par_zoo10->shuffleRSI();//shuffle indices
		par_zoo10->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//GLASS 10
		start = clock();
		clusters_sol = par_glass10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_glass_10 += to_string(i+1) + " Elapsed (Greedy PAR GLASS 10): " + to_string(elapsed) + "(seconds)\n";

		greedy_glass_10 += "\tInfeas: " + to_string(results->Infeasable(par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";
		greedy_glass_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_glass10->atributos, par_glass10->centroides,"GLASS")) + "\n";
		greedy_glass_10 += "\tFitness: " + to_string(results->Fitness(par_glass10->atributos,par_glass10->matriz, clusters_sol, par_glass10->centroides, par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";

		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;

		//cout << i << " GLASS GREEDY 10 ------------------------------" << endl;

		par_glass10->clearClusters(false);//clear the clusters
		par_glass10->shuffleRSI();//shuffle indices
		par_glass10->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//BUPA 10
		start = clock();
		clusters_sol = par_bupa10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_bupa_10 += to_string(i+1) + " Elapsed (Greedy PAR BUPA 10): " + to_string(elapsed) + "(seconds)\n";

		greedy_bupa_10 += "\tInfeas: " + to_string(results->Infeasable(par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";
		greedy_bupa_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa10->atributos, par_bupa10->centroides,"BUPA")) + "\n";
		greedy_bupa_10 += "\tFitness: " + to_string(results->Fitness(par_bupa10->atributos,par_bupa10->matriz, clusters_sol, par_bupa10->centroides, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";

		par_bupa10->clearClusters(false);//clear the clusters
		par_bupa10->shuffleRSI();//shuffle indices
		par_bupa10->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " BUPA GREEDY 10 ------------------------------" << endl;

		//ZOO 20
		start = clock();
		clusters_sol = par_zoo20->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_zoo_20 += to_string(i+1) + " Elapsed (Greedy PAR zoo 20): " + to_string(elapsed) + "(seconds)\n";

		greedy_zoo_20 += "\tInfeas: " + to_string(results->Infeasable(par_zoo20->ML, par_zoo20->CL, par_zoo20->createS())) + "\n";
		greedy_zoo_20 += "\t Error Distance: " + to_string(results->ErrorDistance(clusters_sol, par_zoo20->atributos, par_zoo20->centroides,"ZOO")) + "\n";
		greedy_zoo_20 += "\t Fitness: " + to_string(results->Fitness(par_zoo20->atributos,par_zoo20->matriz, clusters_sol, par_zoo20->centroides, par_zoo20->ML, par_zoo20->CL, par_zoo20->createS())) + "\n";
		//cout << "\tDistance: " << results->Distance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides) << endl;
		par_zoo20->clearClusters(false);//clear the clusters
		par_zoo20->shuffleRSI();//shuffle indices
		par_zoo20->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " ZOO GREEDY 20 ------------------------------" << endl;

		//GLASS 20
		start = clock();
		clusters_sol = par_glass20->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_glass_20 += to_string(i+1) + " Elapsed (Greedy PAR GLASS 20): " + to_string(elapsed) + "(seconds)\n";

		greedy_glass_20 += "\tInfeas: " + to_string(results->Infeasable(par_glass20->ML, par_glass20->CL, par_glass20->createS())) + "\n";
		greedy_glass_20 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_glass20->atributos, par_glass20->centroides,"GLASS")) + "\n";
		greedy_glass_20 += "\tFitness: " + to_string(results->Fitness(par_glass20->atributos,par_glass20->matriz, clusters_sol, par_glass20->centroides, par_glass20->ML, par_glass20->CL, par_glass20->createS())) + "\n";

		par_glass20->clearClusters(false);//clear the clusters
		par_glass20->shuffleRSI();//shuffle indices
		par_glass20->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " GLASS GREEDY 20 ------------------------------" << endl;

		//BUPA 20
		start = clock();
		clusters_sol = par_bupa20->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_bupa_20 += to_string(i+1) + " Elapsed (Greedy PAR BUPA 20): " + to_string(elapsed) + "(seconds)\n";

		greedy_bupa_20 += "\tInfeas: " + to_string(results->Infeasable(par_bupa20->ML, par_bupa20->CL, par_bupa20->createS())) + "\n";
		greedy_bupa_20 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa20->atributos, par_bupa20->centroides,"BUPA")) + "\n";
		greedy_bupa_20 += "\tFitness: " + to_string(results->Fitness(par_bupa20->atributos,par_bupa20->matriz, clusters_sol, par_bupa20->centroides, par_bupa20->ML, par_bupa20->CL, par_bupa20->createS())) + "\n";

		par_bupa20->clearClusters(false);//clear the clusters
		par_bupa20->shuffleRSI();//shuffle indices
		par_bupa20->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " BUPA GREEDY 20 ------------------------------" << endl;

	}
	end_global = clock();
	elapsed = float(end_global - start_global)/CLOCKS_PER_SEC;
	cout << " Elapsed Total (GREEDY): " << elapsed << "(seconds)\n";

	cout << "Calculate BL (ZOO, GLASS, BUPA) Aproximate 7 minutes" << endl;
	start_global = clock();
	for(int i = 0; i < 5; ++i){
		//LOCAL SEARCH
		//ZOO 10

		par_zoo10->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_zoo10->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_zoo_10 += to_string(i+1) + " Elapsed (BL PAR ZOO 10): " + to_string(elapsed) + "(seconds)\n";

		bl_zoo_10 += "\tInfeas: " + to_string(results->Infeasable(par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";
		bl_zoo_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_zoo10->atributos, par_zoo10->centroides,"ZOO")) + "\n";
		bl_zoo_10 += "\tFitness: " + to_string(results->Fitness(par_zoo10->atributos,par_zoo10->matriz, clusters_sol, par_zoo10->centroides, par_zoo10->ML, par_zoo10->CL, par_zoo10->createS())) + "\n";

		par_zoo10->clearClusters(false);//clear the clusters
		par_zoo10->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " ZOO BL 10 ------------------------------" << endl;

		//GLASS 10

		par_glass10->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_glass10->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_glass_10 += to_string(i+1) + " Elapsed (BL PAR GLASS 10): " + to_string(elapsed) + "(seconds)\n";

		bl_glass_10 += "\tInfeas: " + to_string(results->Infeasable(par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";
		bl_glass_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_glass10->atributos, par_glass10->centroides,"GLASS")) + "\n";
		bl_glass_10 += "\tFitness: " + to_string(results->Fitness(par_glass10->atributos,par_glass10->matriz, clusters_sol, par_glass10->centroides, par_glass10->ML, par_glass10->CL, par_glass10->createS())) + "\n";

		par_glass10->clearClusters(false);//clear the clusters
		par_glass10->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " GLASS BL 10 ------------------------------" << endl;

		//BUPA 10
		par_bupa10->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_bupa10->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_bupa_10 += to_string(i+1) + " Elapsed (BL PAR BUPA 10): " + to_string(elapsed) + "(seconds)\n";

		bl_bupa_10 += "\tInfeas: " + to_string(results->Infeasable(par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";
		bl_bupa_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa10->atributos, par_bupa10->centroides,"BUPA")) + "\n";
		bl_bupa_10 += "\tFitness: " + to_string(results->Fitness(par_bupa10->atributos,par_bupa10->matriz, clusters_sol, par_bupa10->centroides, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";

		par_bupa10->clearClusters(false);//clear the clusters
		par_bupa10->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " BUPA BL 10 ------------------------------" << endl;

		//ZOO 20

		par_zoo20->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_zoo20->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_zoo_20 += to_string(i+1) + " Elapsed (BL PAR ZOO 20): " + to_string(elapsed) + "(seconds)\n";

		bl_zoo_20 += "\tInfeas: " + to_string(results->Infeasable(par_zoo20->ML, par_zoo20->CL, par_zoo20->createS())) + "\n";
		bl_zoo_20 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_zoo20->atributos, par_zoo20->centroides,"ZOO")) + "\n";
		bl_zoo_20 += "\tFitness: " + to_string(results->Fitness(par_zoo20->atributos,par_zoo20->matriz, clusters_sol, par_zoo20->centroides, par_zoo20->ML, par_zoo20->CL, par_zoo20->createS())) + "\n";

		par_zoo20->clearClusters(false);//clear the clusters
		par_zoo20->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " ZOO BL 20 ------------------------------" << endl;

		//GLASS 20

		par_glass20->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_glass20->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_glass_20 += to_string(i+1) + " Elapsed (BL PAR GLASS 20): " + to_string(elapsed) + "(seconds)\n";

		bl_glass_20 += "\tInfeas: " + to_string(results->Infeasable(par_glass20->ML, par_glass20->CL, par_glass20->createS())) + "\n";
		bl_glass_20 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_glass20->atributos, par_glass20->centroides,"GLASS")) + "\n";
		bl_glass_20 += "\tFitness: " + to_string(results->Fitness(par_glass20->atributos,par_glass20->matriz, clusters_sol, par_glass20->centroides, par_glass20->ML, par_glass20->CL, par_glass20->createS())) + "\n";

		par_glass20->clearClusters(false);//clear the clusters
		par_glass20->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " GLASS BL 20 ------------------------------" << endl;

		//BUPA 20
		par_bupa20->randomAssign();//assign each node to a cluster randomly

		start = clock();
		clusters_sol = par_bupa20->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_bupa_20 += to_string(i+1) + " Elapsed (BL PAR BUPA 20): " + to_string(elapsed) + "(seconds)\n";

		bl_bupa_20 += "\tInfeas: " + to_string(results->Infeasable(par_bupa20->ML, par_bupa20->CL, par_bupa20->createS())) + "\n";
		bl_bupa_20 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa20->atributos, par_bupa20->centroides,"BUPA")) + "\n";
		bl_bupa_20 += "\tFitness: " + to_string(results->Fitness(par_bupa20->atributos,par_bupa20->matriz, clusters_sol, par_bupa20->centroides, par_bupa20->ML, par_bupa20->CL, par_bupa20->createS())) + "\n";

		par_bupa20->clearClusters(false);//clear the clusters
		par_bupa20->shuffleRSI();
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)
			clusters_sol[i].clear();
		clusters_sol.clear();

		//cout << i << " BUPA BL 20 ------------------------------" << endl;
	}

	end_global = clock();
	elapsed = float(end_global - start_global)/CLOCKS_PER_SEC;
	cout << " Elapsed Total (BL): " << (elapsed/60) << "(minutes)\n" << endl;

	cout << "GREEDY: " << endl;
	cout << endl << "ZOO ********************************************* " << endl;
	cout  << greedy_zoo_10 << endl;
	cout << endl << greedy_zoo_20 << endl;

	cout << "GLASS ********************************************* " << endl;

	cout << endl << greedy_glass_10 << endl;
	cout << endl << greedy_glass_20 << endl;

	cout << "BUPA ********************************************* " << endl;

	cout << endl << greedy_bupa_10 << endl;
	cout << endl << greedy_bupa_20 << endl;

	cout << "LOCAL SEARCH: " << endl;
	cout << endl << "ZOO ********************************************* " << endl;
	cout << bl_zoo_10 << endl;
	cout << endl << bl_zoo_20 << endl;

	cout << "GLASS ********************************************* " << endl;

	cout << endl << bl_glass_10 << endl;
	cout << endl << bl_glass_20 << endl;

	cout << "BUPA ********************************************* " << endl;

	cout << endl << bl_bupa_10 << endl;
	cout << endl << bl_bupa_20 << endl;*/

	/*greedy_bupa_10 = greedy_bupa_20 = bl_bupa_10 = bl_bupa_20 = "";

	par_bupa10->clearClusters(false);//clear the clusters
	par_bupa10->shuffleRSI();//shuffle indices
	par_bupa10->resetCentroides();//reset centroides
	for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
		clusters_sol[i].clear();
	clusters_sol.clear();

	cout << "Calculando Experimento: " << endl;
	//PRUEBA CON UNA SOLUCIÓN OPTIMA GREEDY AL ALGORITMO BL
	for(int i = 0; i<5; ++i){

		//BUPA 10
		start = clock();
		clusters_sol = par_bupa10->algoritmoGreedy();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		greedy_bupa_10 += to_string(i+1) + " Elapsed (Greedy PAR BUPA 10): " + to_string(elapsed) + "(seconds)\n";

		greedy_bupa_10 += "\tInfeas: " + to_string(results->Infeasable(par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";
		greedy_bupa_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa10->atributos, par_bupa10->centroides,"BUPA")) + "\n";
		greedy_bupa_10 += "\tFitness: " + to_string(results->Fitness(par_bupa10->atributos,par_bupa10->matriz, clusters_sol, par_bupa10->centroides, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";


		cout << i << " BUPA GREEDY 10 ------------------------------" << endl;
		start = clock();
		clusters_sol = par_bupa10->algoritmoBL();
		end = clock();
		elapsed = float(end - start)/CLOCKS_PER_SEC;
		bl_bupa_10 += to_string(i+1) + " Elapsed (BL PAR BUPA 10): " + to_string(elapsed) + "(seconds)\n";

		bl_bupa_10 += "\tInfeas: " + to_string(results->Infeasable(par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";
		bl_bupa_10 += "\tError Distance: " + to_string(results->ErrorDistance(clusters_sol, par_bupa10->atributos, par_bupa10->centroides,"BUPA")) + "\n";
		bl_bupa_10 += "\tFitness: " + to_string(results->Fitness(par_bupa10->atributos,par_bupa10->matriz, clusters_sol, par_bupa10->centroides, par_bupa10->ML, par_bupa10->CL, par_bupa10->createS())) + "\n";

		cout << i << " BUPA BL 10 ------------------------------" << endl;

		par_bupa10->clearClusters(false);//clear the clusters
		par_bupa10->shuffleRSI();//shuffle indices
		par_bupa10->resetCentroides();//reset centroides
		for(unsigned int i = 0; i < clusters_sol.size(); ++i)//and clear the clusters solution
			clusters_sol[i].clear();
		clusters_sol.clear();

	}

	cout << "GREEDY BUPA ********************************************* " << endl;

	cout << endl << greedy_bupa_10 << endl;

	cout << "BL BUPA ********************************************* " << endl;

	cout << endl << bl_bupa_10 << endl;*/


	PAR_GM *par_zoo10_gm = new PAR_GM(paths_data[0], paths_const[0], seeds[0]);
	/*PAR_GM *par_zoo20_gm = new PAR_GM("datos/zoo_set.dat", "datos/zoo_set_const_20.const", seed);
	PAR_GM *par_glass10_gm = new PAR_GM("datos/glass_set.dat", "datos/glass_set_const_10.const", seed);
	PAR_GM *par_glass20_gm = new PAR_GM("datos/glass_set.dat", "datos/glass_set_const_20.const", seed);
	PAR_GM *par_bupa10_gm = new PAR_GM("datos/bupa_set.dat", "datos/bupa_set_const_10.const", seed);
	PAR_GM *par_bupa20_gm = new PAR_GM("datos/bupa_set.dat", "datos/bupa_set_const_20.const", seed);
	vector<vector<int>> clusters_sol, clusters_sol_gred, clusters_sol_bl;
	ResultAlgorithms *results = new ResultAlgorithms();

	string greedy_zoo_20 = "", greedy_zoo_10 = "", greedy_glass_10 = "", greedy_glass_20 = "", greedy_bupa_10 = "", greedy_bupa_20 = "";
	string bl_zoo_20 = "", bl_zoo_10 = "", bl_glass_10 = "", bl_glass_20 = "", bl_bupa_10 = "", bl_bupa_20 = "";
	float elapsed;
	clock_t start;
	clock_t end;
	clock_t start_global;
	clock_t end_global;*/

	par_zoo10_gm->randomAssign(50);
	//par_zoo10_gm->printS();
	//par_zoo10_gm->GENETIC(AGG_SF, 0.7);
	par_zoo10_gm->GENETIC(AGG_UN, 0.7, iterations);
	/*par_zoo10_gm->GENETIC(AGE_SF);
	par_zoo10_gm->GENETIC(AGE_UN);*/

	return 0;
}
