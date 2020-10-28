#include "Randoms.cpp"
#include "../include/aco_seq_algorithm.h"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <limits>
#include <climits>
#include <chrono>
#include <ctime>
#include <sstream>
#include <stdlib.h>
#include <fstream>


using namespace std;

Randoms *randoms;

#define ALPHA 1
#define BETA 2
#define Q 50
#define RO 0.5
#define TAUMAX 2

int n_ants = 0;
double best_fitness = 999999.0;
double* best_solution;

double* n_objects_total;
double* n_objects_types;
double* bpp_items;
double* bin_capacity;

double* pheromones;
double *delta_pheromones;
double *step_probabilities;

double *paths;
double *fitnesses;

bool is_palma = false;


// ===  FUNCTION  ======================================================================
//         Name:  run_aco_bpp
//         Description:  EXECUTION ENTRY POINT
// =====================================================================================
double run_aco_bpp(int n_ant, int n_iterations, int problem_id, int isPalma){

	//reset variable
	best_fitness = 999999.0;

	//Execution Time measure
	double average_execution_time = 0.0;
	int iteration = 0;

	//Init structures and variables
	if(isPalma == 0){
		init( n_ant, problem_id, n_iterations, false);
	}else{
		init( n_ant, problem_id, n_iterations, true);
	}

	//START iterations
	while(iteration < n_iterations){
		//packing
		packing();

		//pheromone evaporation
		pheromone_evaporation();

		//Pheromone Deposit
		pheromone_deposit();

		iteration++;
	}

	free(n_objects_types);
	free(n_objects_total);
	free(bin_capacity);
	free(bpp_items);
	free(best_solution);
	free(paths);
	free(fitnesses);

	return best_fitness;
}

// ===  FUNCTION  ======================================================================
//         Name:  INIT()
//         Description:  Create graph, connect cities and generate pheromone matrix
// =====================================================================================
double init(int nAnts, int problem, int iterations, bool palma) {

	is_palma = palma;
	n_ants = nAnts;

	n_objects_types = (double*)malloc(sizeof(double));
	n_objects_total = (double*)malloc(sizeof(double));
	bin_capacity = (double*)malloc(sizeof(double));

	//READ INPUT FILES TO SET MATRIX SIZES
	readBPPFileProperties(problem, n_objects_types, bin_capacity);

	bpp_items = (double*)malloc(2*n_objects_types[0]*sizeof(double)); // 2x -> (weight, quantity)

//	//READ INPUT FILES TO FILL MATRIX SIZES
	readBPPFile(problem, n_objects_types, n_objects_total, bin_capacity, bpp_items);
//
	//alloc other variables
	int bin_capacity_size = bin_capacity[0];
	int n_object_type = n_objects_types[0];
	int n_object_total = n_objects_total[0];
	int pheromone_matrix_size = n_objects_types[0] * n_objects_types[0];

	pheromones = new double[pheromone_matrix_size];
	delta_pheromones = new double[pheromone_matrix_size];
	best_solution = (double*)malloc(n_object_total*sizeof(double));

	paths = (double*)malloc(n_ants* n_object_total * sizeof(double));
	fitnesses = (double*)malloc(n_ants * sizeof(double));

	randoms = new Randoms(15);

	//Create Structures
	initializePheromoneMatrix(n_object_type, pheromones); //Phero OK

	return 0.0;
}

// ===  FUNCTION  ======================================================================
//         Name:  initializePheromoneMatrix
//         Description:  Set Random values to the pheromone matrix
// =====================================================================================
void initializePheromoneMatrix(int n_objects, double* phero){

	double randn = 0.0;

	for(int j = 0;j<n_objects;j++){
		for(int k = 0;k<n_objects;k++){
			if(j!= k){
				randn = randoms -> Uniforme() * TAUMAX;
				phero[(j*n_objects) + k] = randn;
				phero[(k*n_objects) + j] = randn;
			}
			else{
				phero[(j*n_objects) + k] = 0.0;
				phero[(k*n_objects) + j] = 0.0;
			}
		}
	}
}

// ===  FUNCTION  ======================================================================
//         Name:  packing
//         Description:
// =====================================================================================
void packing(){

	double* bpp_items_copy = (double*)malloc(2*((int)n_objects_types[0])*sizeof(double)); // 2x -> (weight, quantity)

	//Iterate over each Ant
	for(int ant_index = 0 ; ant_index < n_ants ; ant_index++){

		//parameters
		int bins_used = 0;
		double n_obj_type = n_objects_types[0];
		double n_obj_total = n_objects_total[0];
		double bin_cap = bin_capacity[0];

		int paths_index = ant_index * n_obj_total;

		//Copy Object and Quantities array
		for(int i = 0 ; i < (int)2*n_obj_type ; i++){
			bpp_items_copy[i] = bpp_items[i];
		}
		
		//Start to compute first bin
		int n_items_in_actual_bin = 0;
		double actual_bin_weight = 0.0;

		double* solution = (double*)malloc(n_obj_total*sizeof(double)); // 2x -> (weight, quantity)
		
		//Used to check if there are still objects that could fit in the actual bin
		int possible_items_to_this_bin = 0;

		//Start first bin -> Get heaviest item available and add to first bin
		int object_index = 0;
		double object_weight = 0.0;
		double new_object_quantity = 0.0;
		double new_object_weight = 0.0;
		
		//Find Heaviest object
		for(int i = 0 ; i < n_obj_type; i++){

			new_object_weight = bpp_items_copy[2*i];
			new_object_quantity = bpp_items_copy[2*i+1];

			if((new_object_quantity > 0) && (new_object_weight > object_weight)){
				object_index = i;
				object_weight = new_object_weight;
			}
		}
		
		//printf("\n 1 - packing Item %i      - Weight %f ", object_index, object_weight);

		//Set first object and update measures
		solution[0] = object_index;
		bpp_items_copy[object_index*2+1]--;
		n_items_in_actual_bin++;
		actual_bin_weight += object_weight;
		bins_used++;

		paths[paths_index] = object_index;

		double* eta = (double*)malloc(n_obj_type*sizeof(double));
		double* tau = (double*)malloc(n_obj_type*sizeof(double));
		double* probs = (double*)malloc(n_obj_type*sizeof(double));

		//Loop to build complete bins
		for (int i = 0; i < n_obj_total-1; i++) {

			double eta_tau_sum = 0.0;

			//Loop to check the possibility of adding other objects
			for (int j = 0; j < n_obj_type; j++) {

				eta[j] = 0.0;
				tau[j] = 0.0;
				probs[j] = 0.0;

				//Get data from the object list
				int weight_object_j = bpp_items_copy[2*j];
				int quantity_object_j = bpp_items_copy[2*j+1];

				//Check if there is still objects available and if the weight suits actual bin
				if((quantity_object_j > 0) && (weight_object_j < (bin_cap-actual_bin_weight))){

					//Calculate the first part of the probability calculation
					if(actual_bin_weight == 0){
						eta[j] = 1;
					}else{
						for(int k = 0 ; k < n_items_in_actual_bin ; k++){
							//last item added to by this ant = index + 1
							//Stay inside last bin using k
							int object_i = solution[i-k];

							eta[j] += pheromones[object_i*(int)n_obj_type+ j];
						}
						eta[j] = eta[j] / n_items_in_actual_bin;
					}

					//Calculate the second part of the probability calculation
					tau[j] = (double) pow(weight_object_j, BETA);

					eta_tau_sum += eta[j] * tau[j];
					possible_items_to_this_bin++;
				}
			}//End checking

			if(possible_items_to_this_bin > 0){
				//Loop to Calculate probabilities based on the values calculated above
				for (int j = 0; j < n_obj_type; j++) {

					probs[j] = (eta[j] * tau[j]) / eta_tau_sum;

					//printf("\n PROBS Object %i, PROB %f", j, d_probs[object_bin_index+j]);

					//Reset Values
					eta[j] = 0.0;
					tau[j] = 0.0;
				}

				//Reset Value
				eta_tau_sum = 0.0;

				//Add new object in a probabilistic manner
				double random = randoms -> Uniforme();
				int object_j = 0;

//				printf("\n Random %f", random);

				double sum = probs[0];
	//			printf("\n sum %f", sum);
				while (sum < random){
					object_j++;
					sum += probs[object_j];
	//				printf("\n sum %f", sum);
				}

				//Add selected object to the list
				solution[i+1] = object_j;
				paths[paths_index+i+1] = object_j;

				//Add weight to actual bin
				double weight_object_j = bpp_items_copy[2*object_j];
				actual_bin_weight += weight_object_j;

				//printf("\n %i - packing Item %i      - Weight %f ",i+2, object_j, weight_object_j);

				//Remove one available item + Increase items in Bin + reset number of possible items to the bin.
				bpp_items_copy[2*object_j+1]--;
				n_items_in_actual_bin++;
				possible_items_to_this_bin = 0;
			}else{

	//			printf("\n\n New BIN ");
				//Start new BIN
				possible_items_to_this_bin = 0;
				actual_bin_weight = 0.0;
				actual_bin_weight = 0;

				//Start first bin -> Get heaviest item available and add to first bin
				int object_index = 0;
				double object_weight = 0.0;
				double object_quantity = 0.0;
				double new_object_weight = 0.0;

				for(int k = 0 ; k < n_obj_type ; k++){

					object_quantity = bpp_items_copy[2*k+1];
					new_object_weight = bpp_items_copy[2*k];

					if((object_quantity > 0) && (new_object_weight > object_weight)){
						object_index = k;
						object_weight = new_object_weight;
					}
				}

				bpp_items_copy[object_index*2+1]--;
				solution[i+1] = object_index;
				paths[paths_index+i+1] = object_index;
				n_items_in_actual_bin++;
				actual_bin_weight += object_weight;

				//printf("\n %i - packing Item %i      - Weight %f ",i+2, object_index, object_weight);

				bins_used++;
			}
		}//end solution

		fitnesses[ant_index] = bins_used;

		if(bins_used < best_fitness){
			best_fitness = bins_used;
//			printf("\n New Best Solution:  %f bins \n", best_fitness);

			for(int l = 0 ; l < n_obj_total; l++){
				best_solution[l] = solution[l];
//				printf(" %f ", best_solution[l]);
			}
		}

		free(solution);
		free(eta);
		free(tau);
		free(probs);
	}//end ants

	free(bpp_items_copy);
}//End PAcking

// ===  FUNCTION  ======================================================================
//         Name:  readBPPFileProperties
//         Description:  read BPP file containing the information about how many product
//			exist in order to create the right structure for the objects.
// =====================================================================================
void pheromone_evaporation(){

	double pheromone_matrix_size = n_objects_types[0] * n_objects_types[0];

	for(int i = 0; i < pheromone_matrix_size ; i++){

		double phero = pheromones[i];
		phero = (1-RO)*phero;
		pheromones[i] = phero;
	}
}//END PHEROMONE EVAPORATION

// ===  FUNCTION  ======================================================================
//         Name:  readBPPFileProperties
//         Description:  read BPP file containing the information about how many product
//			exist in order to create the right structure for the objects.
// =====================================================================================
void pheromone_deposit(){

	double pheromone_matrix_size = n_objects_types[0] * n_objects_types[0];

	int n_obj_total = (int)n_objects_total[0];
	double bin_cap = bin_capacity[0];

	for(int ant_index = 0; ant_index < n_ants; ant_index++){

		double actual_bin_weight = 0.0;
		int actual_bin_n_objects = 0; //final object from bin
		int actual_bin_object_index = 0; //initial object

		for(int i = 0 ; i < n_obj_total ; i++){

		double object_weight = paths[ant_index*n_obj_total+i];

		if(actual_bin_weight + object_weight < bin_cap){
			actual_bin_n_objects++;
			actual_bin_weight+=object_weight;
		}else{
			//update pheromones between items from actual bin index -> n-objects
			for(int j = 0; j<actual_bin_n_objects; j++){
				for(int k = j+1; k<actual_bin_n_objects; k++){

					int object_i = paths[ant_index*n_obj_total+actual_bin_object_index+j];
					int object_j = paths[ant_index*n_obj_total+actual_bin_object_index+k];

					double delta_pheromone =  Q / fitnesses[ant_index];

					pheromones[object_i * (int)n_objects_types[0] + object_j] += delta_pheromone;
					pheromones[object_j * (int)n_objects_types[0] + object_i] += delta_pheromone;
				}
			}

			//Start new bin count
			actual_bin_n_objects = 1;
			actual_bin_weight = object_weight;
			actual_bin_object_index = i;
		}
	}

	}
}//END PHEROMONE EVAPORATION


// ===  HOST FUNCTIONS =================================================================
// =====================================================================================
// ===  FUNCTION  ======================================================================
//         Name:  readBPPFileProperties
//         Description:  read BPP file containing the information about how many product
//			exist in order to create the right structure for the objects.
// =====================================================================================
void readBPPFileProperties(int problem, double* n_objects_type, double* bin_capacity){
//	printf("\n\n Reading BPP File");
	std::ifstream fileReader;

	//Problem Instances
	std::string palma_path = "";

	if(is_palma){
		palma_path = "/home/b/b_mene01/ls_pi-research_aco-bpp/BPP/LowLevelProgram/build/release/";
	}

	//Problem Instances
	std::string f60 = "Falkenauer_t60_00.txt";
	std::string p201 = "201_2500_NR_0.txt";
	std::string p402 = "402_10000_NR_0.txt";
	std::string p600 = "600_20000_NR_0.txt";
	std::string p801 = "801_40000_NR_0.txt";
	std::string p1002 = "1002_80000_NR_0.txt";

	switch(problem){
	case 0:
		fileReader.open(palma_path+f60, std::ifstream::in);
		break;
	case 1:
		fileReader.open(palma_path+p201, std::ifstream::in);
		break;
	case 2:
		fileReader.open(palma_path+p402, std::ifstream::in);
		break;
	case 3:
		fileReader.open(palma_path+p600, std::ifstream::in);
		break;
	case 4:
		fileReader.open(palma_path+p801, std::ifstream::in);
		break;
	case 5:
		fileReader.open(palma_path+p1002, std::ifstream::in);
		break;
	default:
		break;
	}


	if (fileReader.is_open()) {

		fileReader >> n_objects_type[0];
		fileReader >> bin_capacity[0];
	}

//	printf("\n %f --- %f", n_objects_type[0],bin_capacity[0]);

	fileReader.close();
}


// ===  FUNCTION  ======================================================================
//         Name:  readBPPFile
//         Description:  read BPP file containing the problem instance values
// =====================================================================================
void readBPPFile(int problem, double* n_objects_type, double* n_objects_total, double* bin_capacity, double* items){

	std::ifstream fileReader;

	std::string palma_path = "";
	if(is_palma){
		palma_path = "/home/b/b_mene01/ls_pi-research_aco-bpp/BPP/LowLevelProgram/build/release/";
	}

	//Problem Instances
	std::string f60 = "Falkenauer_t60_00.txt";
	std::string p201 = "201_2500_NR_0.txt";
	std::string p402 = "402_10000_NR_0.txt";
	std::string p600 = "600_20000_NR_0.txt";
	std::string p801 = "801_40000_NR_0.txt";
	std::string p1002 = "1002_80000_NR_0.txt";

	switch(problem){
	case 0:
		fileReader.open(palma_path+f60, std::ifstream::in);
		break;
	case 1:
		fileReader.open(palma_path+p201, std::ifstream::in);
		break;
	case 2:
		fileReader.open(palma_path+p402, std::ifstream::in);
		break;
	case 3:
		fileReader.open(palma_path+p600, std::ifstream::in);
		break;
	case 4:
		fileReader.open(palma_path+p801, std::ifstream::in);
		break;
	case 5:
		fileReader.open(palma_path+p1002, std::ifstream::in);
		break;
	default:
		break;
	}

	int lines = 0;
	double total = 0.0;

	if (fileReader.is_open()) {

		fileReader >> n_objects_type[0];
		fileReader >> bin_capacity[0];

		while (lines < n_objects_type[0] && !fileReader.eof()) {
			double weight;
			double quantity;

			fileReader >> weight;
			fileReader >> quantity;

			items[lines*2] = weight;
			items[lines*2+1] = quantity;

			total+=quantity;

			lines++;
		}
	}
	 else{
		 printf("\n File not opened");
	}

	n_objects_total[0] = total;

//	printf("\n Object Types: %f" , n_objects_type[0]);
//	printf("\n Object Total: %f" , n_objects_total[0]);

	fileReader.close();
}

// ===  FUNCTION  ======================================================================
// PRINT FUNCTIONS:
// =====================================================================================

//void printconnections() {
//	cout << " connections: " << endl;
//	cout << "  | ";
//	for (int i = 0; i < n_cities; i++) {
//		cout << i << " ";
//	}
//	cout << endl << "- | ";
//	for (int i = 0; i < n_cities; i++) {
//		cout << "- ";
//	}
//	cout << endl;
//	int count = 0;
//	for (int i = 0; i < n_cities; i++) {
//		cout << i << " | ";
//		for (int j = 0; j < n_cities; j++) {
//			if (i == j) {
//				cout << "x ";
//			} else {
//				cout << connections[i * n_cities + j] << " ";
//			}
//			if (connections[i * n_cities + j] == 1) {
//				count++;
//			}
//		}
//		cout << endl;
//	}
//	cout << endl;
//	cout << "Number of connections: " << count << endl << endl;
//}
//void printRESULTS() {
//	best_length += distance(best_route[n_cities - 1], initial_city);
//	cout << " BEST ROUTE:" << endl;
//	for (int i = 0; i < n_cities; i++) {
//		cout << best_route[i] << " ";
//	}
//	cout << endl << "length: " << best_length << endl;
//
////	cout << endl << " IDEAL ROUTE:" << endl;
////	cout << "0 7 6 2 4 5 1 3" << endl;
////	cout << "length: 127.509" << endl;
//}
//
//void printPOSITIONS() {
//	printf("POSITIONS");
//
//	for (int i = 0; i < n_cities; i++) {
//		printf("\n City : %i", i);
//		printf(" X : %f", cities[i * 3 + 1]);
//		printf(" Y : %f", cities[i * 3 + 2]);
//	}
//}
//
//void printpheromones() {
//	cout << " pheromones: " << endl;
//	cout << "  | ";
//	for (int i = 0; i < n_cities; i++) {
//		printf("%5d   ", i);
//	}
//	cout << endl << "- | ";
//	for (int i = 0; i < n_cities; i++) {
//		cout << "--------";
//	}
//	cout << endl;
//	for (int i = 0; i < n_cities; i++) {
//		cout << i << " | ";
//		for (int j = 0; j < n_cities; j++) {
//			if (i == j) {
//				printf("%5s   ", "x");
//				continue;
//			}
//			if (exists(i, j)) {
//				printf("%7.3f ", pheromones[i * n_cities + j]);
//			} else {
//				if (pheromones[i * n_cities + j] == 0.0) {
//					printf("%5.0f   ", pheromones[i * n_cities + j]);
//				} else {
//					printf("%7.3f ", pheromones[i * n_cities + j]);
//				}
//			}
//		}
//		cout << endl;
//	}
//	cout << endl;
//}
