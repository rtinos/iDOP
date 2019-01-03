/******************************************************************************\
*								 Definitions							 *
\******************************************************************************/
#include <iostream>
using namespace std; 

/* Data structures */
typedef int allele; 									// data structure alele
typedef struct {
			allele *chromosome;							// chromosome
			double fitness;								// fitness
} individual;
typedef struct {			
			individual *ind;
			double sum_fitness;
			double mean_fitness;
			double max_fitness;
			int best_individual;	
			int popsize;	
} population;

// Global variables
extern population popold, popnew, popmem;				// populations
extern double *file_best_fitness, *time_run;			// data to be stored
extern double *file_best_fitness_gen, *global_optimum;	// data to be stored
extern long int gen;									// generation
extern long int max_gen;								// maximum number of generations (only to save the data)
extern int popsize;										// size of the population 
extern int tau;											// change frequency: number of generations to change the landscape
extern double rho;										// change severity
extern int change_type;									// change type: 0:no change, 1:type 1.1 , 2:type 1.2 , 3:type 1.3 , 4:type 2.1 , 5:type 2.2 ,  6:type 2.3 ,7:type 3, 8:mixed
extern int change_cycle_nc;								// define the cycle of change: define the number of changes util stationary landscape
extern int stationary_landscape;						// 1 if the landscape is stationary and 0 otherwise
extern double rimig_rate;								// percentage of the population replaced by random immigrants
extern double mem_rate;									// percentage of the population replaced by immigrants from population memory
extern int n_runs_max ;									// runs of the GA
extern double p_cross;									// crossover rate
extern int tour_pool_size;								// size of the pool for  tournament selection 
extern int tournament_size;								// tournament size
extern int lcrom;										// size of the chromosome 
extern int *knap_w, *knap_p, knap_C;					// knapsack problem: weights, profits, capacity
extern double knap_rho_pen2;							// knapsack problem: penalty parameter

// Declaration of the functions
// statistics.cpp
void statistics( population *pop );
// util_functions.cpp
int *aloc_vectori(int lines);
double *aloc_vectord(int lines);
individual *aloc_vectorind(int lines);
int **aloc_matrixi(int lines , int collums);
double **aloc_matrixd(int lines , int collums);
void desaloc_matrixi(int **Matrix , int lines);
void desaloc_matrixd(double **Matrix , int lines);
int random_int(int L_range, int H_range);
double random_dou(void);
void rand_perm(int *inp, int *out, int size);
void XOR(int *v1, int *v2, int *v3, int l);
int binvec2dec(int *x, int l);
// selection.cpp
int selection( population *pop );
// transformation.cpp
void crossover(allele *parent1, allele *parent2, allele *offspring1, allele *offspring2);
void mutation (allele *offspring, double p_mut);
// fitness.cpp
void initFitKnap( void );
void endFitKnap(void);
double compFitKnap( allele *ind );
double dynProgKnap( void );
// file_man.cpp
void file_output();
//immigrants.cpp
int insertRandomImmig( int i_begin );
void updatePopMem(void);
int insertMemImmig(  int i_begin );

