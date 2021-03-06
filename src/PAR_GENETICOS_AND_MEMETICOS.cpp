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

	//create landa
	landa = createLanda();
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

int PAR_GM::betterFitness(const vector<int> &chromosom, int gen, int &iterations){
	int betterCluster = chromosom[gen], i=0;
	float actual_f, better_f=fitness(chromosom);
	vector<int> copy_chromosom = chromosom;

	//if the value of gen is 0 increase i
	if(i == chromosom[gen]) ++i;

	while(i<k){

		//change cluster
		copy_chromosom[gen] = i;
		//calculate fitness
		actual_f = fitness(copy_chromosom);

		++iterations;

		//and save the cluster with the best fitness
		if(actual_f < better_f){
			better_f = actual_f;
			betterCluster = i;
		}
		++i;
		//if the cluster is equal to actual cluster increment
		if(i == chromosom[gen])
			++i;
	}

	return betterCluster;
}

vector<int> PAR_GM::BL_SOFT(const vector<int> &chromosom, int max_fails, int &iteraciones, int stop){
	vector<int> mejor_sol=chromosom;
	int fails = 0;
	unsigned int i=0;
	bool mejora = true;

	//shuffle the inices
	shuffleRSI();

	//while the solution improve or not exceed the maximum nuber of failures and hasn't covered all solution
	while((mejora or fails < max_fails) and i< RSI.size() and iteraciones < stop){
		mejora = false;

		mejor_sol[RSI[i]] = betterFitness(mejor_sol, RSI[i], iteraciones);

		//if the cluster has changed, there is improvement
		if(chromosom[RSI[i]] != mejor_sol[RSI[i]])
			mejora = true;
		//otherwise it increases the number of faults
		else
			++fails;

		++i;
		++iteraciones;
	}
	return mejor_sol;
}

vector<int> PAR_GM::AM(float probability, int generations, int stop, bool mejores){
	vector<int> mejor_solucion;
	vector<vector<int>> vector_poblacion = vector_solutions;
	vector<float> fitness_poblacion;
	vector<float> fitness_hijos;
	vector<vector<int>> vector_hijos;
	int mejor_padre = -1, peor_hijo= -1, cont_gen=0, n=0, i=0, mejor_hijo=-1;//, iterations=0;
	float mejor_f=999, f_actualp, peor_f = -999, f_actualh;
	vector<int> indices_poblacion;
	vector<vector<int>> vector_mejores;

	//calculate the fitness of actual population
	for(unsigned int e=0; e < vector_poblacion.size(); ++e){

		fitness_poblacion.push_back(fitness(vector_poblacion[e]));

		//and add indices of all chromosomas if not choose 10%*n_better
		if(!mejores)
			indices_poblacion.push_back(e);

		++i;
	}

	//selec the parents
	vector<vector<int>> vector_padres = selectionOperator(vector_poblacion, vector_poblacion.size(),fitness_poblacion);


	//iterate 100.000 evaluations
	while(i<stop){

		//++iterations;
		//choose this cross, because have better results than fixed segment cross
		vector_hijos = uniformMutation(uniformCross(vector_padres, 0.7));
		//vector_hijos = uniformMutation(fixedSegmentCross(vector_padres, probability));

		//if reach 10 generations calculate BL SOFT
		if(cont_gen == 10 && mejores){

			if(indices_poblacion.size()>0)
				indices_poblacion.clear();

			//calculate the number of chromosomes to calculate BL SOFT
			n = vector_hijos.size()*probability;

			for(int e = 0; e<n && i<stop; ++e){
				for(unsigned int j = 0; j<vector_hijos.size(); ++j){
					if(e>0){

						if(fitness_hijos[j] < mejor_f && find(indices_poblacion.begin(), indices_poblacion.end(), j) == indices_poblacion.end()){
							mejor_f = f_actualh;
							mejor_hijo = j;
						}
					}else{

						f_actualh = fitness(vector_hijos[j]);
						fitness_hijos.push_back(f_actualh);
						if(f_actualh < mejor_f){
							mejor_f = f_actualh;
							mejor_hijo = j;
						}
					}
				}

				indices_poblacion.push_back(mejor_hijo);
				mejor_hijo = -1;
				mejor_f = 999;
			}

			//and calculate BL SOFT
			for(unsigned int e=0; e<indices_poblacion.size() && i<stop; ++e)
				vector_hijos[indices_poblacion[e]] = BL_SOFT(vector_hijos[indices_poblacion[e]], 0.1*vector_hijos[indices_poblacion[e]].size(), i, stop);

			//reset the count of generations
			cont_gen = 0;
		}else if(cont_gen == 10){
			//shuffle the new population
			random_shuffle(indices_poblacion.begin(), indices_poblacion.end());//barajo el vector

			//calculate the number of chromosomes to calculate BL SOFT
			n = indices_poblacion.size()*probability;

			//and calculate BL SOFT
			for(int e=0; e<n && i<stop; ++e)
				vector_hijos[indices_poblacion[e]] = BL_SOFT(vector_hijos[indices_poblacion[e]], 0.1*vector_hijos[indices_poblacion[e]].size(), i, stop);

			//reset the count of generations
			cont_gen = 0;
		}

		//choose the best parent and worst descendent
		for(unsigned int e=0; e < vector_poblacion.size(); ++e){

			f_actualp = fitness_poblacion[e];
			f_actualh = fitness(vector_hijos[e]);

			//update fitness of the new population
			fitness_poblacion[e] = f_actualh;

			if(f_actualp < mejor_f){
				mejor_padre = e;
				mejor_f = f_actualp;
			}

			if(f_actualh > peor_f){
				peor_hijo = e;
				peor_f = f_actualh;
			}
			++i;
		}

		//replace the worst descendent for the best parent
		vector_hijos[peor_hijo] = vector_poblacion[mejor_padre];
		fitness_poblacion[peor_hijo] = mejor_f;

		//update vector with new parents
		vector_poblacion = vector_hijos;

		//choose the new parents
		vector_padres = selectionOperator(vector_poblacion, vector_poblacion.size(), fitness_poblacion);

		//reset
		mejor_f = 999;
		mejor_padre = -1;
		peor_f = -999;
		peor_hijo = -1;

		//count number of generations
		++cont_gen;
	}

	mejor_f = 999;
	//choose the best solution
	for(unsigned int e=0; e < vector_poblacion.size(); ++e){

		if(fitness_poblacion[e] < mejor_f){
			mejor_f = fitness_poblacion[e];
			mejor_solucion = vector_poblacion[e];
		}
	}

	return mejor_solucion;
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
	int parent1=-1, parent2=-1, i=0;

	//calculate fitness and the worst 2 parents
	for(unsigned int e=0; e < vector_solutions.size(); ++e){

		f_actual = fitness(vector_poblacion[e]);
		vector_fitness.push_back(f_actual);

		if(f_actual > f_peor1){
			parent2 = parent1;
			f_peor2 = f_peor1;
			parent1 = e;
			f_peor1 = f_actual;
		}else if(f_actual > f_peor2){
			parent2 = e;
			f_peor2 = f_actual;
		}
		++i;
	}

	while(i<stop){

		// select the best parents of each tournament
		vector_padres = selectionOperator(vector_poblacion, 2, vector_fitness);

		//choose the type of cross and calculate the mutation
		if(cruce == AGE_UN)
			vector_hijos = uniformMutation(uniformCross(vector_padres, probability));
		else
			vector_hijos = uniformMutation(fixedSegmentCross(vector_padres, probability));

		//change the worts parents by new descendents
		vector_poblacion[parent1] = vector_hijos[0];
		vector_poblacion[parent2] = vector_hijos[1];

		vector_fitness[parent1] = fitness(vector_poblacion[parent1]);
		vector_fitness[parent2] = fitness(vector_poblacion[parent2]);
		i += 2;

		f_peor1 = f_peor2 = -999;

		//choose the 2 worst parent
		for(unsigned int e=0; e < vector_fitness.size(); ++e){

			if(vector_fitness[e] > f_peor1){
				parent2 = parent1;
				f_peor2 = f_peor1;
				parent1 = e;
				f_peor1 = vector_fitness[e];
			}else if(vector_fitness[e] > f_peor2){
				parent2 = e;
				f_peor2 = vector_fitness[e];
			}
		}
	}

	return vector_poblacion;
}

vector<vector <int>> PAR_GM::AGG(TIPE_CROSS cruce, float probability, int stop){

	vector<vector<int>> vector_poblacion = vector_solutions;
	vector<vector<int>> vector_hijos;
	vector<float> fitness_poblacion;
	int mejor_padre = -1, peor_hijo= -1, i=0;
	float mejor_f=999, f_actualp, peor_f = -999, f_actualh;

	//calculate the fitness of actual population
	for(unsigned int e=0; e < vector_poblacion.size(); ++e){

		fitness_poblacion.push_back(fitness(vector_poblacion[e]));
		++i;
	}

	// select the best parents of each tournament
	vector<vector<int>> vector_padres = selectionOperator(vector_poblacion, vector_solutions.size(), fitness_poblacion);

	while(i<stop){

		//choose the type of cross and calculate the mutation
		if(cruce == AGG_UN)
			vector_hijos = uniformMutation(uniformCross(vector_padres, probability));
		else
			vector_hijos = uniformMutation(fixedSegmentCross(vector_padres, probability));

		//choose the best parent and worst descendent
		for(unsigned int e=0; e < vector_poblacion.size(); ++e){

			f_actualp = fitness_poblacion[e];
			f_actualh = fitness(vector_hijos[e]);

			fitness_poblacion[e] = f_actualh;

			if(f_actualp < mejor_f){
				mejor_padre = e;
				mejor_f = f_actualp;
			}

			if(f_actualh > peor_f){
				peor_hijo = e;
				peor_f = f_actualh;
			}
			++i;
		}

		//and that parent isn't replaced
		vector_hijos[peor_hijo] = vector_poblacion[mejor_padre];
		fitness_poblacion[peor_hijo] = mejor_f;

		//update vector with new parents
		vector_poblacion = vector_hijos;

		//choose the new parents
		vector_padres = selectionOperator(vector_poblacion, vector_solutions.size(), fitness_poblacion);

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
vector<vector<int>> PAR_GM::selectionOperator(const vector<vector<int>> &actual, int tourney, const vector<float> &fitness_p){
	//choose 2 solutions
	int indv1= -1, indv2= -1, indv11=-1, indv12=-1;
	//save the best solution
	vector<vector<int>> padres;

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

		//and save the best of the 2
		if(fitness_p[indv1] < fitness_p[indv2])
			padres.push_back(actual[indv1]);
		else
			padres.push_back(actual[indv2]);
		if(fitness_p[indv11] < fitness_p[indv12])
			padres.push_back(actual[indv11]);
		else
			padres.push_back(actual[indv12]);
	}

	return padres;
}

//uniform crossover operator
vector<vector<int>> PAR_GM::uniformCross(const vector<vector<int>> &padres,float probability){
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

		descendiente1.clear();
		descendiente2.clear();
		descendiente1.resize(RSI.size(), -1);
		descendiente2.resize(RSI.size(), -1);
	}
	return descendientes;
}

//fixed segment crossover operator
vector<vector<int>> PAR_GM::fixedSegmentCross(const vector<vector<int>> &padres, float probability){
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
vector<vector<int>> PAR_GM::uniformMutation(const vector<vector<int>> &padres){
	vector<vector<int>> descendientes = padres;
	vector<int> clust;
	int gen, new_value, crom=-1;

	//calculate the probability of mutations
	float probability =0.1/RSI.size();
	//calculate the number of genes to mutate
	int n_genes = RSI.size()*padres.size() * probability;

	//mutates n_genes
	for(int i = 0; i<n_genes; ++i){

		//choose random chromosome
		crom = rand() % descendientes.size();

		do{
			//choose random gen
			gen = rand() % RSI.size();
			//create new value 0 to k-1
			new_value = rand() % k;

			//calculates the gene vector of that cluster
			clust = findInCluster(descendientes[crom],descendientes[crom][gen]);

		//while the number of genes is less 2 or the new value is equal to the actual gene
		}while(clust.size() < 2 || new_value == descendientes[crom][gen]);

		//change to the new value
		descendientes[crom][gen] = new_value;
	}

	//return the new population
	return descendientes;
}

//Update the distance
vector<vector<float>> PAR_GM::updateDistance(const vector<int> &nodes){
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

//same the infeasibility(int clust, int actual) except it receives the solution
int PAR_GM::infeasibility(const vector<int> &S_cop){
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
float PAR_GM::distanciaEuclidea(const vector<float> &nod1, const vector<float> &nod2){
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
	vector<int> S;

	size_sol = n;
	if(vector_solutions.size()>0){
		for(unsigned int i = 0; i < vector_solutions.size(); ++i)
			vector_solutions[i].clear();
		vector_solutions.clear();
	}
	vector_solutions.resize(n);

	S.resize(RSI.size());

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

float PAR_GM::fitness(const vector<int> &solution){

	float distance = 0, intra_cluster = 0;
	vector<int> elements;
	vector<vector<float>> centroides;
	centroides = updateDistance(solution);

	//walk through each cluster
	for(int i = 0; i< k; ++i){

		elements = findInCluster(solution,i);
		//sumatorry(euclidean distance of all nodes in the cluster)
		for(vector<int>::iterator it = elements.begin(); it != elements.end(); ++it){
			distance += distanciaEuclidea(atributos[(*it)],centroides[i]);
		}
		//mean intra-cluster distance
		distance = distance / elements.size();
		intra_cluster += distance;
		distance = 0;
	}
	float general_deviation = intra_cluster/k;

	float f = general_deviation + infeasibility(solution) * landa;
	return f;
}

//calculates wich of the 2 individuals is the best
int PAR_GM::betterFitness(const vector<vector<int>> &padres, int indv1, int indv2){

	float f1 = 0, f2 = 0;

	//calculate the fitness indv1
	f1 = fitness(padres[indv1]);

	//calculate the fitness indv2
	f2 = fitness(padres[indv2]);

	//if the new solution is better than current solution
	if(f1 < f2)
		return indv1;//and return the actual iteration

	return indv2;
}

//calculate the general deviation
float PAR_GM::generalDeviation(const vector<int> &s_cop){
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
vector<int> PAR_GM::findInCluster(const vector<int> &s_cop, int clust){
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
float PAR_GM::ErrorDistance(const vector<int> &solution, string type_data_file){
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
