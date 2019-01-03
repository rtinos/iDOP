/******************************************************************************\
*								 Immigrants									   *
\******************************************************************************/
#include "defs.h"
#include <cstdlib>


/******************************************************************************\
*			Inserting random immigrants		   								   *
\******************************************************************************/
int insertRandomImmig( int i_begin ){
	int num_ind, size_imig_pool, gene, i_end;
	
	size_imig_pool=(int) (rimig_rate*popsize);
	i_end=size_imig_pool+i_begin-1;
	if (i_end>popsize-1)
		i_end=popsize-1;	
	for (num_ind=i_begin;num_ind<=i_end;num_ind++)
		for (gene=0;gene<lcrom ;gene++) 
     		popold.ind[num_ind].chromosome[gene] = random_int(0,1);
     
	 return(i_end+1);			
}

/******************************************************************************\
*			Updating memory population 	   						   				*
\******************************************************************************/
void updatePopMem(void){
	int gene, i;
				
	if (popmem.popsize<popsize){
		i=popmem.popsize;
		popmem.popsize=popmem.popsize+1;
	}
	else {
		i=random_int(0,popmem.popsize-1);
		while (i==popmem.best_individual)
			i=random_int(0,popmem.popsize-1);
	}		
	for (gene=0;gene<lcrom ;gene++) 
		popmem.ind[i].chromosome[gene] = popold.ind[popold.best_individual].chromosome[gene];
	popmem.ind[i].fitness = popold.ind[popold.best_individual].fitness;
	if (popmem.popsize==1 || popmem.ind[i].fitness>popmem.ind[popmem.best_individual].fitness)
		popmem.best_individual=i;						
				
}


/******************************************************************************\
*			Inserting immigrants from memory population    					   *
\******************************************************************************/
int insertMemImmig(  int i_begin ){
	int gene, i, size_mem_pool, num_ind, i_end;
	
	if (popmem.popsize>mem_rate*popsize)			
		size_mem_pool=(int) (mem_rate*popsize);
	else
		size_mem_pool=popmem.popsize;	
	i_end=size_mem_pool+i_begin-1;
	if (i_end>popsize-1)
		i_end=popsize-1;	
	
	for (num_ind=i_begin;num_ind<=i_end;num_ind++)
		if (num_ind==i_begin){						
			for (gene=0;gene<lcrom ;gene++) 
	     			popold.ind[num_ind].chromosome[gene] = popmem.ind[popmem.best_individual].chromosome[gene];
	    }
	    else{
	     	i=random_int(0,popmem.popsize-1);
	     	for (gene=0;gene<lcrom ;gene++) 
	     		popold.ind[num_ind].chromosome[gene] = popmem.ind[i].chromosome[gene];	     				
		}
					
		return(i_end+1);	
}
