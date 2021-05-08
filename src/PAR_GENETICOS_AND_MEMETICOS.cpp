/*
 * lector_ficheros.cpp
 *
 *  Created on: 23 april 2021
 *      Author: Ruben Girela Castellón
 */
#include <iostream>
#include <armadillo>// para la matriz
#include <fstream> //lectura de ficheros
#include <stdlib.h>
#include <string.h>
#include <math.h> // para calculos matematicos (sqrt, pow (x^y), ...)
#include <algorithm>// para barajar el vector
#include <ctime>// libreria de time para el shuffle aleatorio de vectores
#include <stack>// use stack

#include "../include/PAR_GENETICOS_AND_MEMETICOS.h"
#include "../include/random.h"
using namespace std;
using namespace arma;

PAR_GM::PAR_GM(string fichero_set, string fichero_set_const, int semilla){

	//creo el vector de atributos y la matriz de restricciones
	lectura(fichero_set, fichero_set_const);

	//en función del numero de atributos que tenga el numero de clusters sera 16 o 7
	switch(atributos[0].size()){
	case 5:
		k=16;
		break;
	default:
		k=7;
		break;
	}

	//insert all indices of atributos
	for(unsigned int i=0; i<atributos.size(); ++i){
		RSI.push_back(i);
	}

	//despues barajo los indices de los atributos
	srand(unsigned (semilla));//genero una semilla fija
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR_GM::shuffleRSI(){
	random_shuffle(RSI.begin(), RSI.end());//barajo el vector
}

void PAR_GM::lectura(string fichero_set, string fichero_set_const){
	//leo el fichero *_set.dat
	ifstream read(fichero_set);

	//si no lo encuentra o no puede habrirlo termina el programa
	if (!read){
		cout<<"No puedo abrir "<<fichero_set<<endl;
		exit(1);
	}

	//lee el contenido del archivo
	string cad="", sub_cad="";
	vector<float> fila;
	int size=0;//guarda el tamaño de la matriz
	while(!read.eof()) {
		++size;
		read >> cad;//leo por filas
		//recorro la fila caracter por caracter
		for(unsigned int i=0; i<cad.size(); ++i){
			if(cad.at(i)!=','){//si no encuentra una coma guarda el numero
				sub_cad+=cad.at(i);
			}else{//si encuentra una coma es que ya tiene el numero entero
				fila.push_back(stof(sub_cad));//lo añado al vector auxiliar
				sub_cad="";//y reseteo la sub cadena
			}
		}
		//añado el ultimo numero que no tiene delimitador a la derecha
		fila.push_back(stof(sub_cad));
		sub_cad="";//reseteo la subcadena

		//añado el vector de distancias de la fila leida
		atributos.push_back(fila);
		fila.clear();//limpio el vector
	}
	read.close();
	//borro el ultimo vector, ya que es la del salto de linea de la ultima linea, y contiene los mismos valores que la anterior
	atributos.pop_back();

	--size;//le decremento 1 porque cuenta el salto de linea del vacio

	matriz.set_size(size,size);//redimensiono la matriz
	size_mat = size;

	read.open(fichero_set_const);//abro el archivo que contiene la matriz

	//si no lo encuentra o no puede habrirlo termina el programa
	if (!read){
		cout<<"No puedo abrir "<<fichero_set_const<<endl;
		exit(1);
	}

	//lee el contenido del archivo *_set_const
	cad="";
	string num="";
	int col=0, row=0, col2=1;
	while(!read.eof()) {
		read >> cad;

		if(row<size){//esto lo hago porque lee el salto de linea de la ultima linea
			//recorro la linea
			for(unsigned int i=0; i<cad.size(); i+=2){
				//guardo el numero
				num = cad.at(i);
				if(num=="-"){//si es negativo guardo el - y el numero
					++i;
					num += cad.at(i);
				}
				//solo guardo la diagonal
				if(col>=col2){
					//guardo el valor en la matriz
					matriz(row,col)=stod(num);

					if(col != row and matriz(row,col) != 0){
						if(matriz(row,col) >0)
							ML.push_back(pair<int,int>(row,col));
						if(matriz(row,col) < 0)
							CL.push_back(pair<int,int>(row,col));
					}
				}

				++col;
			}
			//incremento de la fila e inicialización de la columna
			++row;
			col=0;
			col2=row+1;
		}

	}
	read.close();
}

//reset Centroides
/*void PAR_GM::resetCentroides(){
	//vuelvo a generar aleatoriamente los centroides inicialmente distintos de dimensión n
	for(int i=0; i<k; ++i){
		centroides[i].clear();
		for(unsigned int e=0; e <atributos[0].size(); ++e)
			centroides[i].push_back(Rand());
	}
}*/

//print the elements of each cluster
void PAR_GM::printS(int s){
	vector<int> elements;
	int count = 0, total = 0;
	if(s == 0){

		for(int j=0; j< size_sol; ++j){
			cout << "SOLUCIÓN: " << j+1 << endl;
			for(int i = 0; i<k; ++i){
				elements = findInCluster(vector_solutions[j],i);
				cout << i << ": [ ";
				for(auto e = elements.begin(); e != elements.end(); ++e){
					cout << *e << ", ";
					++count;
				}
				cout << " ] n = " << count << endl;
				total += count;
				count = 0;
			}
			cout << "total elements = " << total << endl;
			total = 0;
		}
	}else{
		for(int i = 0; i<k; ++i){
			elements = findInCluster(vector_solutions[s],i);
			cout << i << ": [ ";
			for(auto e = elements.begin(); e != elements.end(); ++e){
				cout << *e << ", ";
				++count;
			}
			cout << " ] n = " << count << endl;
			total += count;
			count = 0;
		}
		cout << "total elements = " << total << endl;
	}
}

//print the elements of the solution
void PAR_GM::printSolution(vector<int> s){
	vector<int> elements;
	int count = 0, total = 0;

	for(int i = 0; i<k; ++i){
		elements = findInCluster(s,i);
		cout << i << ": [ ";
		for(auto e = elements.begin(); e != elements.end(); ++e){
			cout << *e << ", ";
			++count;
		}
		cout << " ] n = " << count << endl;
		total += count;
		count = 0;
	}
	cout << "total elements = " << total << endl;
}

//imprime los atributos de cada nodo
void PAR_GM::printDistanciasEuclideas(){
	for(vector<vector<float>>::iterator it=atributos.begin(); it != atributos.end(); ++it){
		for(vector<float>::iterator it2=(*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << endl;
	}
}
//imprime los centroides de cada cluster
void PAR_GM::printCentroides(vector<vector<float>> centroides){
	int i=0;
	for(vector<vector<float>>::iterator it=centroides.begin(); it != centroides.end(); ++it){
		cout << i << " [ ";
		for(vector<float>::iterator it2=(*it).begin(); it2 != (*it).end(); ++it2){
			if(it2+1 != (*it).end())
				cout << (*it2) << ", ";
			else
				cout << (*it2);
		}
		cout << " ]" << endl;
		++i;
	}
}

//print the indexes of the attributes
void PAR_GM::printRSI(){
	for(vector<int>::iterator it=RSI.begin(); it != RSI.end(); ++it){
		if(it+1 != RSI.end())
			cout << (*it) << ", ";
		else
			cout << (*it) << endl;
	}
}

vector<int> PAR_GM::GENETIC(TIPE_CROSS tipo, float probability, int stop){

	vector<vector<int>> mejores;
	vector<int> mejor_solucion;
	float f_actual, mejor_f = 999;

	//choose the tipe of crossover operator and genetic algorithm
	if(tipo == AGG_UN || tipo == AGG_SF){
		mejores = AGG(tipo, probability, stop);
	}else
		mejores = AGE(tipo, probability, stop);

	for(unsigned int i=0; i < mejores.size(); ++i){

		f_actual = fitness(mejores[i]);

		if(f_actual < mejor_f){
			mejor_f = f_actual;
			mejor_solucion = mejores[i];
		}
	}

	return mejor_solucion;
}

vector<vector <int>> PAR_GM::AGE(TIPE_CROSS cruce, float probability, int stop){

	vector<vector<int>> vector_poblacion = vector_solutions;
	vector<vector<int>> vector_padres;
	vector<vector<int>> vector_hijos;
	vector<float> vector_fitness;
	float f_actual, f_peor1=-999, f_peor2=-999;
	int parent1=-1, parent2=-1, iterations = 0;

	for(unsigned int i=0; i < vector_solutions.size(); ++i){
		f_actual = fitness(vector_poblacion[i]);
		vector_fitness.push_back(f_actual);
		if(f_actual > f_peor1){
			parent2 = parent1;
			f_peor2 = f_peor1;
			parent1 = i;
			f_peor1 = f_actual;
		}else if(f_actual > f_peor2){
			parent2 = i;
			f_peor2 = f_actual;
		}
	}

	for(int i=0; i<stop; i+=2){

		vector_padres = selectionOperator(vector_poblacion, 2);

		if(cruce == AGE_UN)
			vector_hijos = uniformMutation(uniformCross(vector_padres, probability));
		else
			vector_hijos = uniformMutation(fixedSegmentCross(vector_padres, probability));

		vector_poblacion[parent1] = vector_hijos[0];
		vector_poblacion[parent2] = vector_hijos[1];

		vector_fitness[parent1] = fitness(vector_poblacion[parent1]);
		vector_fitness[parent2] = fitness(vector_poblacion[parent2]);

		f_peor1 = f_peor2 = -999;

		//choose the 2 worst parent
		for(unsigned int i=0; i < vector_fitness.size(); ++i){

			if(vector_fitness[i] > f_peor1){
				parent2 = parent1;
				f_peor2 = f_peor1;
				parent1 = i;
				f_peor1 = f_actual;
			}else if(f_actual > f_peor2){
				parent2 = i;
				f_peor2 = f_actual;
			}
		}
		++iterations;
		//cout << iterations << endl;
	}

	//cout << "NUMBER OF ITERATIONS: " << iterations << endl;

	return vector_poblacion;
}

vector<vector <int>> PAR_GM::AGG(TIPE_CROSS cruce, float probability, int stop){

	vector<vector<int>> vector_poblacion = vector_solutions;
	vector<vector<int>> vector_padres = selectionOperator(vector_poblacion, vector_solutions.size());
	vector<vector<int>> vector_hijos;
	int mejor_padre = -1, peor_hijo= -1;
	float mejor_f=999, f_actualp, peor_f = -999, f_actualh;

	stop = stop/vector_solutions.size();
	//cout << stop << endl;

	for(int i=0; i<stop; ++i){
		//cout << "size: " << vector_padres.size() << endl;
		//cout << "iteracion: " << i << endl;

		//choose the type of cross and calculate the mutation
		if(cruce == AGG_UN)
			vector_hijos = uniformMutation(uniformCross(vector_padres, probability));//uniformMutation(uniformCross(vector_padres));
		else
			vector_hijos = uniformMutation(fixedSegmentCross(vector_padres, probability));

		//cout << "iteracion despues: " << i << endl;
		//cout << "size: " << vector_padres.size() << endl;

		/*for(auto p:vector_padres){
			for(auto g:p){
				cout << g << ", ";
			}
			cout << endl;
		}*/
		/*for(auto p:vector_hijos){
			for(auto g:p){
				cout << g << ", ";
			}
			cout << endl;
		}*/

		//cout << "poblacion: " << vector_poblacion.size() << endl;
		//cout << "hjos: " << vector_hijos.size() << endl;

		//choose the best parent
		for(unsigned int e=0; e < vector_poblacion.size(); ++e){

			//cout << "antes fitness: size = " << vector_padres[e].size() << ", e = " << e << endl;
			//cout << "size hijos: " << vector_hijos[e].size() << endl;
			/*for(int j = 0; j<vector_padres[e].size(); ++j){
				cout << vector_hijos[e][j] << ", ";
			}
			cout << endl;*/
			f_actualp = fitness(vector_poblacion[e]);
			f_actualh = fitness(vector_hijos[e]);

			//cout << "despues fitness" << endl;

			if(f_actualp < mejor_f){
				mejor_padre = e;
				mejor_f = f_actualp;
			}

			if(f_actualh > peor_f){
				peor_hijo = e;
				peor_f = f_actualh;
			}

		}
		//cout << "aqui?, hijo: " << peor_hijo << ", padre: " << mejor_padre << endl;
		//and that parent isn't replaced
		vector_hijos[peor_hijo] = vector_poblacion[mejor_padre];

		//cout << "MEJOR PADRE: " << mejor_f << endl;

		//update vector with new parents
		vector_poblacion = vector_hijos;

		//choose the new parents
		vector_padres = selectionOperator(vector_poblacion, vector_solutions.size());

		//reset
		mejor_f = 999;
		mejor_padre = -1;
		peor_f = -999;
		peor_hijo = -1;
	}

	//and return it
	return vector_poblacion;
}
//select the best new set of solutions
vector<vector<int>> PAR_GM::selectionOperator(vector<vector<int>> actual, int tourney){
	//choose 2 solutions
	int indv1= -1, indv2= -1, indv11=-1, indv12=-1;
	//save the best solution
	vector<vector<int>> padres;
	//create landa
	landa = createLanda();

	//generate n torney depend of tipe AGG or AGE
	for(int i=0; i < tourney; i+=2){
		//randomly select 2 diferents individuals
		do{
			indv1 = rand() % actual.size();
			indv2 = rand() % actual.size();
		}while(indv1 == indv2);
		do{
			indv11 = rand() % actual.size();
			indv12 = rand() % actual.size();
		}while(indv11 == indv12);

		//cout << indv1 << " vs " << indv2 << endl;

		//cout << "THE BEST: " << betterFitness(actual,indv1,indv2) << endl;
		//and save the best of the 2
		padres.push_back(actual[betterFitness(actual,indv1,indv2)]);
		padres.push_back(actual[betterFitness(actual,indv11,indv12)]);
	}

	return padres;
}

//uniform crossover operator
vector<vector<int>> PAR_GM::uniformCross(vector<vector<int>> padres,float probability){
	vector<vector<int>> descendientes = padres;
	vector<int> RSI_CROSS = RSI;
	vector<int> descendiente1;
	vector<int> descendiente2;

	//create 2 new descendents
	int number_cross = (probability*padres.size())/2;
	int parent1=-1, parent2=-1;

	/*cout << (probability*padres.size())/2 << endl;
	cout << number_cross << endl;*/

	descendiente1.resize(RSI.size(), -1);
	descendiente2.resize(RSI.size(), -1);

	for(int i=0; i< number_cross; ++i){

		//choose random 2 parents
		do{
			parent1 = rand() % descendientes.size();
			parent2 = rand() % descendientes.size();
		}while(parent1 == parent2);


		//generates n/2 different random indices different from genes one parent and the rest of the other parent
		random_shuffle(RSI_CROSS.begin(), RSI_CROSS.end());//barajo el vector

		int size1 = RSI_CROSS.size()/2;
		int size2 = RSI_CROSS.size() - size1;
		int size_max = size1;

		if(size2 > size1)
			size_max = size2;

		//cout << size1 << " vs " << size2 << " = " << size_max << endl;

		//create the son
		for(int e=0; e<size_max; ++e){

			//n/2 genes first parent
			if(e < size1){
				descendiente1[RSI_CROSS[e]] = descendientes[parent1][RSI_CROSS[e]];
				descendiente2[RSI_CROSS[e]] = descendientes[parent2][RSI_CROSS[e]];
			}
			//and rest of the second parent
			//cout << e+size1 << endl;
			descendiente1[RSI_CROSS[e+size1]] = descendientes[parent2][RSI_CROSS[e+size1]];
			descendiente2[RSI_CROSS[e+size1]] = descendientes[parent1][RSI_CROSS[e+size1]];

		}
		//add the children
		descendientes[parent1] = descendiente1;
		descendientes[parent2] = descendiente2;


		/*cout << "HIJO: " << i+1 << endl;
		for(int j = 0; j<k; ++j){
			vector<int> elements = findInCluster(descendiente1,j);
			cout << j << ": [ ";
			for(auto e : elements){
				cout << e << ", ";
			}
			cout << " ] n = " << elements.size() << endl;
		}
		cout << "HIJO: " << i+2 << endl;
		for(int j = 0; j<k; ++j){
			vector<int> elements = findInCluster(descendiente2,j);
			cout << j << ": [ ";
			for(auto e : elements){
				cout << e << ", ";
			}
			cout << " ] n = " << elements.size() << endl;
		}*/

		descendiente1.clear();
		descendiente2.clear();
		descendiente1.resize(RSI.size(), -1);
		descendiente2.resize(RSI.size(), -1);
	}
	return descendientes;
}

//fixed segment crossover operator
vector<vector<int>> PAR_GM::fixedSegmentCross(vector<vector<int>> padres, float probability){
	vector<vector<int>> descendientes = padres;
	//random start segment and size of segment
	unsigned int start_seg1 = 0, start_seg2 = 0, size_seg1 = 0, end_seg1=0, size_seg2 = 0, end_seg2=0;
	//size of other elements, number of elements to first and second parent
	int size_other1 = 0, size_other2=0, n_parent11 = 0, n_parent12 = 0, n_parent21 = 0, n_parent22 = 0, choose_parent=-1;
	bool not_seg1 = true, not_seg2 = true;
	//create 2 new descendents
	vector<int> descendiente1;
	vector<int> descendiente2;

	//calculate the number of pairs to cross
	int number_cross = (probability*padres.size())/2;
	int parent1=-1, parent2=-1;

	descendiente1.resize(RSI.size(), -1);
	descendiente2.resize(RSI.size(), -1);

	for(int i=0; i< number_cross; ++i){

		//choose random 2 parents
		do{
			parent1 = rand() % descendientes.size();
			parent2 = rand() % descendientes.size();
		}while(parent1 == parent2);

		//generate the start segment
		start_seg1 = rand() % RSI.size();
		start_seg2 = rand() % RSI.size();
		//the size of segment
		size_seg1 = 1 + rand() % (RSI.size() - 1);
		size_seg2 = 1 + rand() % (RSI.size() - 1);
		//and the end of the segment
		end_seg1 = (start_seg1 + size_seg1)%RSI.size() - 1;
		end_seg2 = (start_seg2 + size_seg2)%RSI.size() - 1;
		//size of other descents
		size_other1 = RSI.size() - size_seg1;
		size_other2 = RSI.size() - size_seg2;
		//and divide that size between the 2 parents
		n_parent11 = size_other1 /2;
		n_parent12 = size_other1 - n_parent11;

		n_parent21 = size_other2 /2;
		n_parent22 = size_other2 - n_parent21;

		/*cout << "hijo1 --------------" << endl;
		cout << "START1: " << start_seg1 << endl;
		cout << "END1: " << end_seg1 << endl;
		cout << "SIZE1: " << size_seg1 << endl;
		cout << "SIZE REST1: " << size_other1 << endl;
		cout << "N PARENT11: " << n_parent11 << endl;
		cout << "N PARENT12: " << n_parent12 << endl;
		cout << "TOTAL: " << n_parent11+n_parent12+size_seg1 << endl;
		cout << "hijo2 --------------" << endl;
		cout << "START2: " << start_seg2 << endl;
		cout << "END2: " << end_seg2 << endl;
		cout << "SIZE2: " << size_seg2 << endl;
		cout << "SIZE REST2: " << size_other2 << endl;
		cout << "N PARENT21: " << n_parent21 << endl;
		cout << "N PARENT22: " << n_parent22 << endl;
		cout << "TOTAL: " << n_parent21+n_parent22+size_seg2 << endl;*/


		//go through all elements
		for(unsigned int e = 0; e< RSI.size(); ++e){

			//if element is within range add element of the fisrt parent
			//if the range is 0...k
			if(start_seg1 < end_seg1){

				if(e>= start_seg1 && e<=end_seg1){

					descendiente1[e] = descendientes[parent1][e];
					not_seg1 = false;
				}
			//if the range is k..n and 0..j
			}else{

				if(e>= start_seg1){

					descendiente1[e] = descendientes[parent1][e];
					not_seg1 = false;
				}else if(e<=end_seg1){

					descendiente1[e] = descendientes[parent1][e];
					not_seg1 = false;
				}
			}

			//if element is within range add element of the second parent
			//if the range is 0...k
			if(start_seg2 < end_seg2){

				if(e>= start_seg2 && e<=end_seg2){

					descendiente2[e] = descendientes[parent2][e];
					not_seg2 = false;
				}
			//if the range is k..n and 0..j
			}else{

				if(e>= start_seg2){

					descendiente2[e] = descendientes[parent2][e];
					not_seg2 = false;
				}else if(e<=end_seg2){

					descendiente2[e] = descendientes[parent2][e];
					not_seg2 = false;
				}
			}

			if(not_seg1){

				//choose the parent randomly
				if(n_parent11>0 and n_parent12>0)
					choose_parent = rand() % 2;

				//and add gen of parent in descendent
				//parent 1
				//if you have choose the first parent and still can add genes from that first parent or the second parent  has already selected all the genes
				if((choose_parent == 0 && n_parent11 >0) || n_parent12<=0){
					descendiente1[e] = descendientes[parent1][e];
					--n_parent11;
				//parent 2
				//if you have choose the second parent and still can add genes from that second parent or the first parent  has already selected all the genes
				}else if((choose_parent == 1 && n_parent12 >0) || n_parent11<=0){
					descendiente1[e] = descendientes[parent2][e];
					--n_parent12;
				}
			}
			if(not_seg2){

				//choose the parent randomly
				if(n_parent21>0 and n_parent22>0)
					choose_parent = rand() % 2;

				//and add gen of parent in descendent
				if((choose_parent == 0 && n_parent21 >0) || n_parent22<=0){
					descendiente2[e] = descendientes[parent1][e];
					--n_parent21;
				}else if((choose_parent == 1 && n_parent22 >0) || n_parent21<=0){
					descendiente2[e] = descendientes[parent2][e];
					--n_parent22;
				}
			}
			not_seg1 = true;
			not_seg2 = true;

		}

		/*int total = 0;
		cout << "HIJO 1: " << i+1 << endl;
		for(int j = 0; j<k; ++j){
			vector<int> elements = findInCluster(descendiente1,j);
			cout << j << ": [ ";
			for(auto e : elements){
				cout << e << ", ";
			}
			cout << " ] n = " << elements.size() << endl;
			total += elements.size();
		}
		cout << "total = " << total << endl;
		total = 0;
		cout << "HIJO 2: " << i+2 << endl;
		for(int j = 0; j<k; ++j){
			vector<int> elements = findInCluster(descendiente2,j);
			cout << j << ": [ ";
			for(auto e : elements){
				cout << e << ", ";
			}
			cout << " ] n = " << elements.size() << endl;
			total += elements.size();
		}
		cout << "total = " << total << endl;*/
		//replace the parents with their children
		descendientes[parent1]=descendiente1;
		descendientes[parent2]=descendiente2;

		//and clear children
		descendiente1.clear();
		descendiente1.resize(RSI.size(), -1);
		descendiente2.clear();
		descendiente2.resize(RSI.size(), -1);
	}


	return descendientes;
}

//uniform mutation operator
vector<vector<int>> PAR_GM::uniformMutation(vector<vector<int>> padres){
	vector<vector<int>> descendientes = padres;
	vector<int> clust;
	int gen, new_value, crom=-1;

	//calculate the probability of mutations
	float probability =0.1/RSI.size();
	//calculate the number of genes to mutate
	int n_genes = RSI.size()*padres.size() * probability;

	/*cout << "genes: " << RSI.size() << ", cromosomas: " << padres.size() << endl;
	cout << probability << ", " << n_genes << endl;*/

	//mutates n_genes
	for(int i = 0; i<n_genes; ++i){

		//cout << "GEN: " << i << endl;

		//choose random chromosome
		crom = rand() % descendientes.size();

		do{
			//choose random gen
			gen = rand() % RSI.size();
			//create new value 0 to k-1
			new_value = rand() % k;

			//cout << "GEN: " << gen << ", Value: " << new_value << endl;

			//calculates the gene vector of that cluster
			clust = findInCluster(descendientes[crom],descendientes[crom][gen]);

		//while the number of genes is less 2 or the new value is equal to the actual gene
		}while(clust.size() < 2 || new_value == descendientes[crom][gen]);

		//cout << "(ORIGINAL) GEN: " << gen << ", Value: " << descendientes[crom][gen] << endl;

		//change to the new value
		descendientes[crom][gen] = new_value;
		//cout << "(NEW) GEN: " << gen << ", Value: " << descendientes[crom][gen] << endl;
	}
	//return the new population
	return descendientes;
}

//Update the distance
vector<vector<float>> PAR_GM::updateDistance(vector<int> nodes){
	//save the actual distance
	vector<float> distance(atributos[nodes[0]].size(),0);
	vector<vector<float>> centroides;
	centroides.resize(k);

	int size = 0;
	for(int j = 0; j<k; ++j){
		for(unsigned int e=0; e<nodes.size(); ++e){
			if(nodes[e] == j){
				++size;
				//calculate average distance
				for(unsigned int i=0; i < atributos[e].size(); ++i){
					//sumatorry
					distance[i] += atributos[e][i];
				}
			}
		}

		//average
		for(unsigned int i = 0; i < distance.size(); ++i){
			distance[i] = distance[i]/nodes.size();
		}
		centroides[j] = distance;
	}


	return centroides;
}

//calculate the closest and least restriction to cluster
int PAR_GM::minRestrictionsDistance(int actual, bool first){
	int cluster=-1;
	/*float min_distance=999;//save the minimum distance and less restriction
	int less_restriction=999;
	float actual_distance=0;//save the actual distance
	int actual_restriction=0;//save the actual number of restrictions

	//go through all clusters
	for(unsigned int i=0; i < centroides.size(); ++i){

		//calculate the Euclidea distance with the current cluster
		actual_distance = distanciaEuclidea(atributos[actual],centroides[i]);

		//if it isn't the first node to enter
		if(!first)
			//calculates the number of constraints it violates
			actual_restriction = infeasibility(i, actual);

		//if the current cluster is less than the minimum saved, update the cluster and distance
		if(actual_restriction <=less_restriction){

			if((actual_restriction <less_restriction) || (actual_distance < min_distance && actual_restriction <=less_restriction)){

				min_distance = actual_distance;
				less_restriction = actual_restriction;
				cluster = i;
			}
		}
	}*/
	return cluster;
}

//calculate infeasibility when assigning an atribute to each cluster and return the minimum
int PAR_GM::infeasibility(int clust, int actual){
	// number of restrictions, matrix column and row and not empty cluster indexes
	int rest=0;
	//walk through each constrain ML
	for (unsigned int i = 0; i<ML.size(); ++i) {
		if(i < ML.size()){
			//if it find the current node in the constrain
			if(ML[i].first == actual){
				//and the second node has an assigned cluster and isn't in the same cluster
				if(S[ML[i].second] != -1 && S[ML[i].second] != clust){
					++rest;//increases the number of restrictions violated
				}
			}else if(ML[i].second == actual){
				if(S[ML[i].first] != -1 && S[ML[i].first] != clust){
					++rest;
				}
			}
		}
	}
	//walk through each constrain CL
	for (unsigned int i = 0; i<CL.size(); ++i) {
		//if it find the current node in the constrain
		if(CL[i].first == actual){
			//and the second node has an assigned cluster and is in the same cluster
			if(S[CL[i].second] != -1 && S[CL[i].second] == clust){
				++rest;//increases the number of restrictions violated
			}
		}else if(CL[i].second == actual){
			if(S[CL[i].first] != -1 && S[CL[i].first] == clust){
				++rest;
			}
		}
	}

	return rest;// return the number of restrictions it violates
}

//same the infeasibility(int clust, int actual) except it receives the solution
int PAR_GM::infeasibility(vector<int> S_cop){
	int restrictions = 0;
	//CALCULATE the max size
	unsigned int max = CL.size();
	if(ML.size()>max) max = ML.size();

	//walk through all restrictions
	for(unsigned int i = 0; i< max; ++i){
		if(i<CL.size()){
			//check if it violate the CL constraint
			if(S_cop[CL[i].first] == S_cop[CL[i].second])
				++restrictions;
		}

		if(i<ML.size()){
			//check if it violate the ML constraint
			if(S_cop[ML[i].first] != S_cop[ML[i].second])
				++restrictions;
		}
	}

	return restrictions;
}

//calcula la distancia euclidea entre 2 nodos
float PAR_GM::distanciaEuclidea(vector<float> nod1, vector<float> nod2){
	//la formula es: sqrt(sumatoria((a_i - b_i)²))
	float suma=0;

	//sumatoria((a_i - b_i)²)
	for(unsigned int i=0; i<nod1.size();++i){
		suma += pow((nod1[i]-nod2[i]),2);
	}
	//sqrt(sumatoria)
	return sqrt(suma);
}

//random assignment of each node with a cluster
void PAR_GM::randomAssign(int n){
	bool not_null = true;
	//bool find = false;
	stack<int> clusters_null;
	vector<int> clusters_asign;

	size_sol = n;
	vector_solutions.resize(n);

	S.resize(RSI.size());
	//cout << S.size() << endl;

	//initialize vector
	for(int i = 0; i<n; ++i){
		vector_solutions[i].resize(RSI.size());
	}

	for(int j=0; j<n; ++j){

		//go through all nodes
		for(unsigned int i = 0; i < RSI.size(); ++i){
			//if it has traversed half the nodes
			if(i == RSI.size()/2){
				//check that no cluster is empty
				for(int e = 0; e < k; ++e){
					/*for(unsigned int l=0; l<RSI.size()/2 && !find; ++l){
						if(vector_solutions[j][(int)RSI[l]] == e)
							find=true;
					}*/
					//if the cluster is empty
					/*if(!find){//add the cluster to the stack
						not_null = false;
						clusters_null.push(e);
					}else
						find = false;*/

					//find cluster empty
					if(find(clusters_asign.begin(), clusters_asign.end(), e) == clusters_asign.end()){

						clusters_null.push(e);
						not_null = false;
					}
				}
			}


			//if all clusters aren't empty or haven't yet traveled half of the nodes
			if(not_null){
				vector_solutions[j][RSI[i]] = rand() % k;

				if(find(clusters_asign.begin(), clusters_asign.end(), vector_solutions[j][RSI[i]]) == clusters_asign.end()){
					clusters_asign.push_back(vector_solutions[j][RSI[i]]);
				}

			}else{//else it asigns to the cluster that is empty
				vector_solutions[j][RSI[i]] = clusters_null.top();
				clusters_null.pop();
			}
			//if all the clusters aren't empty, assign the rest of the nodes randomly
			if(clusters_null.size() == 0)
				not_null = true;

		}
	}
}

float PAR_GM::fitness(vector<int> solution){
	float f = generalDeviation(solution) + infeasibility(solution) * landa;
	return f;
}

//calculates wich of the 2 individuals is the best
int PAR_GM::betterFitness(vector<vector<int>> padres, int indv1, int indv2){

	float f1 = 0, f2 = 0;

	//calculate the fitness indv1
	f1 = fitness(padres[indv1]);

	//calculate the fitness indv2
	f2 = fitness(padres[indv2]);

	//cout << "Fitness: " << f1 << " vs " << f2 << endl;

	//if the new solution is better than current solution
	if(f1 < f2)
		return indv1;//and return the actual iteration

	return indv2;
}

//calculate the general deviation
float PAR_GM::generalDeviation(vector<int> s_cop){
	float distance = 0, intra_cluster = 0;
	vector<int> elements;
	vector<vector<float>> centroides;
	centroides = updateDistance(s_cop);

		//walk through each cluster
		for(int i = 0; i< k; ++i){

			elements = findInCluster(s_cop,i);
			//sumatorry(euclidean distance of all nodes in the cluster)
			for(vector<int>::iterator it = elements.begin(); it != elements.end(); ++it){
				distance += distanciaEuclidea(atributos[(*it)],centroides[i]);
			}
			//mean intra-cluster distance
			distance = distance / elements.size();
			intra_cluster += distance;
			distance = 0;
		}

		return intra_cluster/k;
}

//find all elements of the cluster
vector<int> PAR_GM::findInCluster(vector<int> s_cop, int clust){
	vector<int> elements;
	for(unsigned int i=0; i<s_cop.size(); ++i)
		if(clust == s_cop[i])
			elements.push_back(i);

	return elements;
}

//calculate Landa
float PAR_GM::createLanda(){
	float lan = 0, actual_distance=0;

	//calculate the maximum distance
	for(unsigned int i = 0; i < atributos.size(); ++i){
		for(unsigned int e = i+1; e < atributos.size(); ++e){
			actual_distance = distanciaEuclidea(atributos[i],atributos[e]);
			if(actual_distance > lan){
				lan = actual_distance;
			}
		}
	}
	//count the number of restrictions total
	uvec rest = (find(matriz == 1 or matriz == -1));

	//max distance / total number of problem restrictions
	lan = lan/rest.size();

	return lan;
}

//calculate the distance error
float PAR_GM::ErrorDistance(vector<int> solution, string type_data_file){
	float error_distance = 0, optimal_distance = 0;
	//calculate the distance
	error_distance = generalDeviation(solution);

	//choose the optimal distance
	if(type_data_file.compare("ZOO") == 0)
		optimal_distance = 0.904799856;
	else if(type_data_file.compare("GLASS") == 0)
		optimal_distance = 0.364290282;
	else
		optimal_distance = 0.220423749;

	//calculate the distance error
	error_distance = abs(error_distance-optimal_distance);

	return error_distance;
}
