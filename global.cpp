#include "defs.h"

/******************************************************************************\
*				  	Global Variables								 		   *
\******************************************************************************/
// GA
population popold, popnew, popmem;				// populations
int maxgen;										// maximum number of generations 
double *file_best_fitness, *time_run;			// data to be stored
double *file_best_fitness_gen, *global_optimum;	// data to be stored
long int gen;									// generation
long int max_gen=10000;							// maximum number of generations (only to save the data)
int lcrom;										// size of the cromosome 
int popsize = 100;								// size of the population (>1)
int n_runs_max = 50;							// runs of the GA
double p_cross=0.6;								// crossover rate
int tournament_size=3;							// tournament size
double rimig_rate;								// percentage of the population replaced by random immigrants
double mem_rate;								// percentage of the population replaced by immigrants from population memory
// DOP
int tau;										// change frequency: number of generations to change the landscape
double rho;										// change severity
int change_type;								// change type: 0:no change, 1:type 1.1 , 2:type 1.2 , 3:type 1.3 , 4:type 2.1 , 5:type 2.2 ,  6:type 2.3 ,7:type 3, 8:mixed
int change_cycle_nc=1;							// define the cycle of change: define the number of changes util stationary landscape
int stationary_landscape;						// 1 if the landscape is stationary and 0 otherwise
// Knapsack fitness function
int *knap_w, *knap_p, knap_C;					// knapsack problem: weights, profits, capacity
double knap_rho_pen2;	
