/******************************************************************************\
*								 Knapsack Problem							 *
\******************************************************************************/
#include "defs.h"
#include <cstdlib>
#include <cmath>


/******************************************************************************\
*			Dynamic Programming (see Wikipedia)		   *
\******************************************************************************/
double dynProgKnap( void ){
	int i, j, **M, temp;
	double opt_fit;
	
	M = aloc_matrixi(lcrom+1,knap_C+1);
	
 	for (j=0;j<=knap_C;j++)
    	M[0][j]=0;
 	for (i=1;i<=lcrom;i++){ 	
 		for (j=0;j<=knap_C;j++){
 			if ( knap_w[i-1]>j )
 				M[i][j]=M[i-1][j];
 			else {
 				temp=M[i-1][j-knap_w[i-1]] + knap_p[i-1] ;
 				if (M[i-1][j]<temp)
 					M[i][j]=temp;
				else
					 M[i][j]=M[i-1][j];	
 			}
 		}
	 }


	opt_fit=(double) M[lcrom][knap_C];
	
	desaloc_matrixi(M,lcrom+1);
	
	return ( opt_fit );
	
}


/******************************************************************************\
*			Setting the parameters for the Knapsack Fitness Function		   *
\******************************************************************************/
void initFitKnap( void ){
	int gene;
	int sum_w=0.0;
	
	// Memory allocation
	knap_p = aloc_vectori(lcrom);
	knap_w = aloc_vectori(lcrom);
	
	knap_rho_pen2 = (double) knap_p[0]/knap_w[0];
	for (gene=0;gene<lcrom;gene++) {
		knap_p[gene] =  random_int(40,100);		// random int in the interval [40,100]
		knap_w[gene] =  random_int(5,20);		// random int in the interval [5,20]
		sum_w = sum_w + knap_w[gene];
		if ( ((double) knap_p[gene]/knap_w[gene]) > knap_rho_pen2 )
			knap_rho_pen2= (double) knap_p[gene]/knap_w[gene];
	}
	knap_C = (int) (0.5 * sum_w); 	
}

/******************************************************************************\
*		  Knapsack Fitness Function Parameters: desallocate memory 		       *
\******************************************************************************/
void endFitKnap(void){

	delete [] knap_p;
	delete [] knap_w;
	
}

/******************************************************************************\
*			Fitness for the Knapsack problem with Pen2 , Han & Kim, 2000	           *
\******************************************************************************/
double compFitKnap( allele *ind ){
		int gene;
		double sum_p=0.0, sum_w=0.0, pen2;

		// Sum of weights
		for (gene=0;gene<lcrom;gene++) 
			sum_w = sum_w + knap_w[gene]*ind[gene];
		// Sum of Profits
		for (gene=0;gene<lcrom;gene++) 
			sum_p = sum_p + knap_p[gene]*ind[gene];
		
		pen2=knap_rho_pen2*(sum_w-knap_C);
		if (pen2>sum_p)
			pen2=sum_p;		// avoiding negative values
			
		if 	(sum_w<=knap_C)
			return (sum_p);
		else if (pen2<sum_p)
			return (sum_p-pen2);
		else
			return (0.0);

}

