#include "defs.h"
#include <cstdlib>


/******************************************************************************\
*								 2-point Crossover	 						 *
\******************************************************************************/
void Point2X(allele *parent1, allele *parent2, allele *offspring1, allele *offspring2 ){
	
		int p1=0, p2=0, aux, gene;

		// defining the crossover points
		while (p1 == p2) {
			p1 =random_int (0,lcrom-1);	// point 1
			p2 =random_int (0,lcrom-1);	// point 2
		}
		if (p1>p2) {
			aux=p1;
			p1=p2;
			p2=aux;
		}
							 
		// generating the offspring
		for (gene=0;gene<lcrom;gene++) {
			if (gene<p1 || gene>=p2) {
				offspring1[gene] = parent1[gene];
				offspring2[gene] = parent2[gene];	
			}
			else{
				offspring1[gene] = parent2[gene];
				offspring2[gene] = parent1[gene];
			}
		}

}


/******************************************************************************\
*								 Uniform Crossover	 						 *
\******************************************************************************/
void UX(allele *parent1, allele *parent2, allele *offspring1, allele *offspring2  ){
		int aux, gene;

		for (gene=0;gene<lcrom;gene++){
			aux=random_int (0,1);				// mask: define if the gene comes from parent 1 (0) or 2 (1)
			if (aux==0){
				offspring1[gene] = parent1[gene];	
				offspring2[gene] = parent2[gene];		
			}
			else{
				offspring1[gene] = parent2[gene];
				offspring2[gene] = parent1[gene];
			}
		}	
				
}


/******************************************************************************\
*								 Crossover	 												*
\******************************************************************************/
void crossover(allele *parent1, allele *parent2, allele *offspring1, allele *offspring2){
	int gene, crossover_type=1; // 1: UX; 2: 2-point crossover

	if ( random_dou () <= p_cross ){
		if (crossover_type==1)
			UX(parent1, parent2, offspring1, offspring2  );
		else if (crossover_type==2)
			Point2X(parent1, parent2, offspring1, offspring2  );
	}
	else{
		for (gene=0;gene<lcrom;gene++){
			offspring1[gene] = parent1[gene];
			offspring2[gene] = parent2[gene];	
		}		
	}
	
}


/******************************************************************************\
*								 Mutation														   *
\******************************************************************************/
void mutation (allele *offspring, double p_mut){
	int gene;
	
	for (gene=0;gene<lcrom;gene++){
		if ( random_dou () < p_mut ){
			if (offspring[gene]==0)
				offspring[gene]=1;
			else
				offspring[gene]=0;
		}
	}	

}


