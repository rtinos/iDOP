*** A Framework for Inducing Artificial Changes in Pseudo-Boolean Optimization Problems ***

Description: This is the source code for the framework for inducing artificial changes. 
This code was used to generate the results for the experiments with the 0-1 knapsack problem.

Reference:  Tinós, R. & Yang, S. (2019). "A Framework for Inducing Artificial Changes in Optimization Problems", Submitted to Information Sciences.	

Contact: Renato Tinos <rtinos@ffclrp.usp.br>

Running the code: ./iDOP_knap <tau> <rho> <change_type> <rimig_rate> <mem_rate> <N>

Parameters:
<tau>: change frequency, i.e., number of generations to change the fitness landscape
<rho>: change severity (0.0<=rho<=1.0)
<change_type>:change type: 0 - no change, 1 - type 1.1 , 2 - type 1.2 , 3 - type 1.3 , 4 - type 2.1 , 5 - type 2.2 ,  6 - type 2.3 , 7 - type 3, 8 - mixed
<rimig_rate>: controls the percentage of the population replaced by random immigrants (only when immigrants algorithm is used)
<mem_rate>: controls the percentage of the population replaced by immigrants from population memory (only when immigrants algorithm is used)
<N>: dimension (number of items in 0-1 knapsack problem

Example for running the code with 200 itens (change type 1.1, with tau=500, rho=0.001 and standard GA): 

make

./iDOP_knap 500 0.001 1 0.0 0.0 200;


Observation 1: Class dop is given in dop.h 

- Function double dop::transform( int *x ,  int *x_new) : This function is used to transform a solution before its evaluation (and also to generate the fitness difference). 
This function can be used independently from the Search Algorithm.
	
- Function void dop::change( int change_type_p, double rho, double f_range, int *x_template ): This function is used to change the fitness landscape.
			
Observation 2: file global.cpp contains the parameters of the GA (examples: population size and crossover rate).

Observation 3: iDOP_knap generates four files. The files, generated in directory "data/Nx/" (where x is the dimension), are:
The first file (bfg_....dat) contains the best fitness in each generation for run 0. 
The second file (bfi_....dat) contains the fitness of best solutions in each run.
The third file (tim_....dat) contains the time for each run.
The fourth file (opt_....dat) contains the fitness of the global optimum given by the Dynamic Programming algorithm (for comparison).

Observation 4: the running time can be modified by changing the line:
"} while (double( clock() - time_start ) / (double)CLOCKS_PER_SEC < ((double) lcrom/1.0) );" 
in the ga function (idop_knap.cpp).