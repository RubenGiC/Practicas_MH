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

#include "../include/PAR.h"
using namespace std;
using namespace arma;

int main() {
	//mat B(4,5,fill::randu);

	//cout << B << endl;
	//leo los archivos
	//"datos/bupa_set.dat", "datos/bupa_set_const_10.const"
	PAR *par = new PAR("datos/zoo_set.dat", "datos/zoo_set_const_10.const");
	//par->lectura("datos/bupa_set.dat", "datos/bupa_set_const_10.const");
	//par.printDistanciasEuclideas();
	//cout << par.matriz << endl;
	//par->printCentroides();
	//par->printRSI();

	vector<vector<int>> clusters;

	clock_t start = clock();
	clusters = par->algoritmoGreedy();
	clock_t end = clock();
	float elapsed = float(end - start)/CLOCKS_PER_SEC;
	cout << "Elapsed (Greedy): " << elapsed << "(seconds)" << endl;

	cout << "Solution clusters:" << endl;
	int n = 0;
	for(vector<vector<int>>::iterator it = clusters.begin(); it != clusters.end(); ++it){
		cout << n << " [ ";
		for(vector<int>::iterator it2 = (*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << " ]" << endl;
		++n;
	}
	par->landa();
	return 0;
}
