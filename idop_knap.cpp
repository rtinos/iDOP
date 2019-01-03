/* -*- Mode: C; indent-tabs-mode: t; c-basic-offset: 4; tab-width: 4 -*-  */
/*
 * Inducing changes in the 0-1 knapsack problem (iDOP: induced Dynamic Optimization Problem)
 *
 * Copyright (C) 2019  Renato Tinos <rtinos@ffclrp.usp.br>
 * 
 * References:
 *		1) Tinos, R. & Yang, S. (2019). "A Framework for Inducing Artificial Changes in Optimization Problems", Submitted to Information Sciences.
 *
 * iDOP_knap is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * iDOP_knap  is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#include "defs.h"
#include <cstdlib>
#include <cmath>
#include <time.h>
#include "dop.h"				// dop class


/******************************************************************************\
*				  	Print data							 			 *
\******************************************************************************/
void print_data(population *pop ){

	
	cout <<"Generation:"<< gen << endl;
	cout <<"Best individual:"<< pop->best_individual << endl;
	cout <<"Fitness of the best individual:"<< pop->max_fitness << endl;
	cout <<"Mean fitness: "<< pop->mean_fitness << endl;
	
	/*	
	int i, gene;
	for (i=0;i<pop->popsize ;i++) {	
		cout <<"("<< pop->ind[i].fitness<<") " ;
		for (gene=0;gene<lcrom ;gene++) 
			cout << pop->ind[i].chromosome[gene]<<" ";
		cout << endl;
	}*/
}


/******************************************************************************\
*						Setting the parameters for the Fitness Function		   *
\******************************************************************************/
void initFitness( void ){
	initFitKnap();	
}


/******************************************************************************\
*					 Fitness Function: memory desallocation 		   			*
\******************************************************************************/
void endFitness( void ){	
	endFitKnap();		
}


/******************************************************************************\
*								Compute Fitness 							 *
\******************************************************************************/
double compFitness( allele *ind , dop *DOP ){
		int *ind_transf;
		double Fitness, delta_f;

		ind_transf = aloc_vectori(lcrom);
		
		delta_f=DOP->transform(ind, ind_transf);

		Fitness = compFitKnap( ind_transf ) + delta_f;		// compFitKnap: fitness function for the knapsack problem

		
		delete [] ind_transf;

		return Fitness;
}


/******************************************************************************\
*								 Generation of the GA							 *
\******************************************************************************/
void generation( double p_mut, int n_run, dop *DOP){
	int gene, j=0 , parent1, parent2;

	do {	
		// Selection of two parents	
		parent1=selection( &popold );
		parent2=selection( &popold );
		// Crossover
		crossover( popold.ind[parent1].chromosome , popold.ind[parent2].chromosome,  popnew.ind[j].chromosome , popnew.ind[j+1].chromosome );	
		// Mutation
		mutation( popnew.ind[j].chromosome, p_mut );
		mutation( popnew.ind[j+1].chromosome, p_mut );
		// Elitism
		if (j==0){		
			for (gene=0;gene<lcrom;gene++)
				popnew.ind[j].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
			popnew.ind[j].fitness=popold.ind[popold.best_individual].fitness;
		}
		else
			popnew.ind[j].fitness=compFitness( popnew.ind[j].chromosome, DOP );		// fitness computation
		popnew.ind[j+1].fitness=compFitness( popnew.ind[j+1].chromosome, DOP );		// fitness computation
		j = j + 2;	
	} while ( j < popnew.popsize);
	
}
				

/******************************************************************************\
*				  	Initiate Population 					 				 *
\******************************************************************************/
void initiatePop(dop *DOP){
	int gene, num_ind;
	
	// Size of the populations
	popold.popsize=popsize;
	popnew.popsize=popsize;
	popmem.popsize=0;
	
	// Dynamic allocation: populations
	popold.ind = aloc_vectorind(popsize);
	popnew.ind = aloc_vectorind(popsize);
	popmem.ind = aloc_vectorind(popsize);
	
	for (num_ind=0;num_ind<popsize;num_ind++){
		// Dynamic allocation: chromosomes	
		popold.ind[num_ind].chromosome = aloc_vectori(lcrom);
		popnew.ind[num_ind].chromosome = aloc_vectori(lcrom);
		popmem.ind[num_ind].chromosome = aloc_vectori(lcrom);

		// Random Initialization
	 	for (gene=0;gene<lcrom ;gene++) {
     		popold.ind[num_ind].chromosome[gene] = random_int(0,1);
		}
        popold.ind[num_ind].fitness = compFitness( popold.ind[num_ind].chromosome, DOP );								// fitenss computation	
        
	}
	
	statistics( &popold);
	//print_data(&popold);
	
}

/******************************************************************************\
*				   Population: desallocate memory 					 				 *
\******************************************************************************/
void endPop(void){
	int num_ind; 
	
	for (num_ind=0;num_ind<popnew.popsize;num_ind++){
		delete [] popold.ind[num_ind].chromosome;
		delete [] popnew.ind[num_ind].chromosome;
		delete [] popmem.ind[num_ind].chromosome;
	}
	delete [] popold.ind;
	delete [] popnew.ind;
	delete [] popmem.ind;
}

/******************************************************************************\
*				  	Copy Population							 			 *
\******************************************************************************/
void copy_pop( void ){
	int i, gene;
		
	for (i=0;i<popnew.popsize ;i++) {	
		popold.ind[i].fitness=popnew.ind[i].fitness;
		for (gene=0;gene<lcrom ;gene++) 
			popold.ind[i].chromosome[gene]=popnew.ind[i].chromosome[gene];
	}
	
}


/******************************************************************************\
*				  	Run of the GA 			 *
\******************************************************************************/
void ga(int n_run, double p_mut){	
	int gen_init, gene, num_ind, change_type_temp=0, change_type_temp_last, ind_index;
	int change_cycle_init, change_cycle;								// change cycle
	double f_range, max_fit;
	clock_t time_start;
					
	// Initializing

	gen = 0;
	gen_init=gen;
	change_cycle = 1;
	change_cycle_init=change_cycle;
	initFitness( );									// initiate fitness function
	dop *DOP = new dop(lcrom);							// create DOP
	initiatePop(DOP);									// initiate population
	max_fit=popold.max_fitness;
	if (n_run==0)
		file_best_fitness_gen[gen]=popold.max_fitness;
	// Computing f_range (used in DOP generator)
	f_range=popold.max_fitness-popold.mean_fitness;
	if (f_range<0.1)
		f_range=popold.max_fitness;
	// Dynamic Programming
	if (stationary_landscape==1 && rimig_rate==0.0 && mem_rate==0.0)
		global_optimum[n_run]=dynProgKnap();
	time_start=clock();		
	
	// Genetic Algorithm
	do {
		gen = gen + 1; 								// generation index		
		//cout<<"gen: "<<gen<<endl;	
		// Change
		if ( (gen-gen_init)>=tau ){
			change_cycle++;							// change cycle increment	
			gen_init=gen;
			change_type_temp_last=change_type_temp;
			// Modification of the landscape (DOP change)					
			if ( (change_cycle-change_cycle_init)>change_cycle_nc ){
				change_type_temp=0;				// stationary environment
				change_cycle_init=change_cycle;
			}
			else{
				if (change_type==8)
					change_type_temp=random_int(1,7);		// change: random type
				else
					change_type_temp=change_type;	
			}
			//cout<<"Changing the problem....change type: "<<change_type_temp<<", rho= "<<rho<<endl;
			DOP->change(change_type_temp, rho, f_range,  popold.ind[popold.best_individual].chromosome);		// change
			
			// Elitism (case immigrants are inserted)
			if (rimig_rate>0.0 || mem_rate>0.0)
				for (gene=0;gene<lcrom;gene++)
					popold.ind[0].chromosome[gene]=popold.ind[popold.best_individual].chromosome[gene];
			ind_index=1;
			// Inserting random immigrants after change	
			if (rimig_rate>0.0){				
				//cout<<"Inserting random immigrants ...."<<endl;				
				// Generating random solutions	and inserting in the population
				ind_index=insertRandomImmig(ind_index);
			}	
			// Inserting immigrants from population memory after change and updating memory	
			if (mem_rate>0.0){	
				if (change_type_temp_last==0){	
					// Updating immigrants memory population	
					//cout<<"Updating  population memory ...."<<endl;	
					updatePopMem();
				}				
				if (change_type_temp==0 && popmem.popsize>0){	
					// Inserting immigrants from memory
					//cout<<"Inserting immigrants from population memory ...."<<endl;	
					ind_index=insertMemImmig(ind_index);
				}
			}
			// Computing the fitness for the new landscape 		
			for (num_ind=0;num_ind<popold.popsize;num_ind++)
				popold.ind[num_ind].fitness = compFitness( popold.ind[num_ind].chromosome , DOP );	// compute fitness												
			statistics( &popold );			
		}
		
		generation(p_mut, n_run, DOP);
			
		copy_pop();		// popold=popnew
		statistics( &popold );			
		//print_data(&popold);
		// Check maximum fitness if the landscape is stationary
		if (change_type_temp==0)
			if (popold.max_fitness>max_fit)
				max_fit=popold.max_fitness;				
		// save the fitness across the generations (only for the first run)
		if (n_run==0 && gen<max_gen)
			file_best_fitness_gen[gen]=popold.max_fitness;
	//} while (double( clock() - time_start ) / (double)CLOCKS_PER_SEC < ((double) lcrom/6.0) );
	} while (double( clock() - time_start ) / (double)CLOCKS_PER_SEC < ((double) lcrom/1.0) );
	
	// Data to be saved
	time_run[n_run] = double( clock() - time_start ) / (double)CLOCKS_PER_SEC;
	file_best_fitness[n_run] = max_fit;
		
	endPop();		// population: desallocation 
	delete DOP;
	endFitness( );	// fitness function: desallocation
}


/******************************************************************************\
*				  	Main													 *
\******************************************************************************/
int main(int argc , char *argv[])
{
	int n_run,  N_input;
	double p_mut;

	// Arguments
	if( argc < 7) {
		cout<<"Insufficient number of arguments!"<<endl;
		cout<<"Call: iDOP <tau> <rho> <change_type> <rimig_rate> <mem_rate> <N>"<<endl;
		exit(1);
	}
	else{
		tau=atoi(argv[1]);	
		rho=atof(argv[2]);
		change_type=atoi(argv[3]);
		stationary_landscape=0;
		if (tau<=0){		
			tau=2147483647;			//  MAX_INT: tau used for stationary environment
			stationary_landscape=1;
		}
		rimig_rate=atof(argv[4]);
		mem_rate=atof(argv[5]);
		N_input=atoi(argv[6]);
		if (tau<0 || rho<0.0 || change_type<0 || change_type>8 || rimig_rate>1.0 || change_cycle_nc<1 || (rimig_rate+mem_rate)>0.9 || N_input<1 ){
			cout<<"Incorrect arguments!"<<endl;
			cout<<"Call: iDOP  tau>=0  rho>=0.0  0<=change_type<=8  rimig_rate<=1.0 change_cycle_nc<1 mem_rate<=1.0 (rimig_rate+mem_rate)<=0.9 <N> (>0) "<<endl;
			exit(1);
		}
	}	
		
	// Parameters 
	lcrom=N_input; // chromosome lenght
	p_mut=(double) 1.0/lcrom;			// mutation rate                		
	
	// Allocation of vectors for the data to be stored
	file_best_fitness_gen=aloc_vectord(max_gen);
	file_best_fitness=aloc_vectord(n_runs_max);
	time_run=aloc_vectord(n_runs_max);
	global_optimum=aloc_vectord(n_runs_max);

	cout << "\n ***** Genetic Algorithm ****" << endl;
	cout << "N="<<N_input << endl;
	for (n_run=0;n_run<n_runs_max;n_run++) {	
		srand(n_run+1);	// random seed   		
		cout <<"Run:"<< n_run << endl;
		ga(n_run,p_mut);
	}
	file_output();					// save data

	delete [] global_optimum;
	delete [] time_run;
	delete [] file_best_fitness;
	delete [] file_best_fitness_gen;
	
	return 0;
}
