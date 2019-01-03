/******************************************************************************\
*								 Statistics						 *
\******************************************************************************/

#include "defs.h"
#include <cstdlib>


void statistics( population *pop )
{	
	int j;

	pop->sum_fitness = pop->ind[0].fitness; 			// sum of the fitness in the population
	pop->max_fitness = pop->ind[0].fitness;   			// maximum fitness in the population
	pop->best_individual = 0;							// best individual in the population

	for(j=1;j<popsize;j++) {
		pop->sum_fitness = pop->sum_fitness + pop->ind[j].fitness;
		if (pop->ind[j].fitness > pop->max_fitness )	{	
			pop->max_fitness = pop->ind[j].fitness ; 
			pop->best_individual = j;			
		}
	}

	pop->mean_fitness = pop->sum_fitness / popsize; 	// mean fitness in the population



}			
