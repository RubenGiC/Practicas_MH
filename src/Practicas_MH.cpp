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
	//par.print_distancias_euclideas();
	//cout << par.matriz << endl;
	//cout << par->min_distancia() << endl;
	//par->print_centroides();
	//par->print_RSI();
	int p = par->min_distance(par->atributos[0]);

	vector<int> a = {1,2,3};
	/*vector<int> b = {1,5,3};*/
	vector<vector<int>> c;
	vector<vector<int>> d;
	c.push_back(a);
	/*d.push_back(b);*/
	cout << c.size() << endl;
	if(c != d)
		cout << "son iguales" << endl;
	else
		cout << "son distintos" << endl;

	return 0;
}
