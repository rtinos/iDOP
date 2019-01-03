/******************************************************************************\
*				  				 Files Manipulation							 *
\******************************************************************************/

#include "defs.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include<string>
#include<cstring>
#include<fstream>

#define CHAR_LEN 1000

/******************************************************************************\
* 										Save data : end of the simulation						 *
\******************************************************************************/
void file_output(void)
{
	int i, tau_t, imig_flag=0, mem_flag=0;
	FILE *Bfg_file;
	FILE *Bestfit_file;
	FILE *Time_file;
	FILE *Opt_file;
	char *name_p;
	char name[CHAR_LEN];

    name_p = name;

	if (stationary_landscape==1)
		tau_t=0;
	else
		tau_t=tau;
	if (rimig_rate>0.0)
		imig_flag=1;
	if (mem_rate>0.0)
		mem_flag=1;
		
  // Best fitness in each generation for run 0
	sprintf(name,"data/N%d/bfg_%d_%d_%d_%d_%d.dat",lcrom,tau_t,((int) (rho*1000)),imig_flag,mem_flag,change_type);
	if ((Bfg_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfg to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<max_gen;i++) {
		fprintf(Bfg_file,"%.0f ",file_best_fitness_gen[i]);
	}
	fclose(Bfg_file);

    // Best fitness 
	sprintf(name,"data/N%d/bfi_%d_%d_%d_%d_%d.dat",lcrom,tau_t,((int) (rho*1000)),imig_flag,mem_flag,change_type);
	if ((Bestfit_file = fopen(name_p,"w"))==NULL) {
		puts("The file bfi to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<n_runs_max;i++) {
		fprintf(Bestfit_file,"%.0f ",file_best_fitness[i]);
	}
	fclose(Bestfit_file);
	
  	// Time for each run
	sprintf(name,"data/N%d/tim_%d_%d_%d_%d_%d.dat",lcrom,tau_t,((int) (rho*1000)),imig_flag,mem_flag,change_type);
	if ((Time_file = fopen(name_p,"w"))==NULL) {
		puts("The file time to be saved cannot be open \n");
		exit(1);
	}
	for (i=0;i<n_runs_max;i++) {
		fprintf(Time_file,"%.2f ",time_run[i]);
	}
	fclose(Time_file);
	
		// Global optima
	if (stationary_landscape==1 && rimig_rate==0.0 && mem_rate==0.0){	
		sprintf(name,"data/N%d/opt_%d_%d_%d_%d.dat",lcrom,tau_t,((int) (rho*1000)),imig_flag,mem_flag);
		if ((Opt_file = fopen(name_p,"w"))==NULL) {
			puts("The file opt to be saved cannot be open \n");
			exit(1);
		}
		for (i=0;i<n_runs_max;i++) {
			fprintf(Opt_file,"%.0f ",global_optimum[i]);
		}
		fclose(Opt_file);
	}
	
}


