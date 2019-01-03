/*****************************************************************************************************************************************************\
*								 Class: DOP (Continuous Optimization Problems)																									  *	
* Main reference:																																	  *
* Tinos, R. & Yang, S. (2019). "A Framework for Inducing Artificial Changes in Optimization Problems", Submitted to Information Sciences.			  *
* Aditional references:			 				   																									  *
*  Ref. 1: Tinos, R. & Yang, S. (2014). "Analysis of fitness landscape modifications in evolutionary dynamic optimization", Information Sciences, 282.*
*  Ref. 2: TinosR. & Yang, S. (2007). "Continuous dynamic problem generators for evolutionary algorithms", Proc. of IEEE CEC, 236-243.				  *
* 																																					  *
* DOP types (Continuous Optimization Problems):																										  *
* 	DOP Type 1: DOP with rotation 																												      *
* 		DOP Type 1.1 (DOPs with rotation of the candidate solutions defined by candidate solution exchanges of the XOR type)						  *
* 		DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation of x)		  *
* 		DOP Type 1.3 (DOPs with rotation of the candidate solutions defined by decision variable exchanges according to template)					  *
* 	DOP Type 2: Single time-dependent DOP obtained by copy of decision variables																	  *
* 		DOP Type 2.1 (DOPs obtained by copying elements of the decision variables in vector x)														  *
* 		DOP Type 2.2 (DOPs obtained by copying decision variables from a random solution according to a template)									  *
* 		DOP Type 2.3 (DOPs obtained by copying decision variables from the current best solution according to a template)							  *
* 	DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a template													  *	
* 		  																																			  *	
* For Continuous Optimization Problems, the fitness of solution x can be computed in the search algorithm using function:							  *													
* 														  																							  *	
* long double compFitness( long double *x , dop *DOP ){														  										  	  *	
* 		long double *x_transf, Fitness, delta_f;														  											  *	
* 		x_transf = new long double[lcrom];														  													  *	
* 		delta_f=DOP->transform(x, x_transf);														  											      *	
* 		Fitness = staticFitness( x_transf ) + delta_f;																  							      *	
* 		delete [] x_transf;																  							      							  *	
* 		return Fitness;																  							      								  *	
* }																					  							      								  *	
* 																					  							      								  *	
\*****************************************************************************************************************************************************/



#include <cmath>
#include <cstdlib>

class dop{
	private:	
		int l;															// size of the chromosome, i.e, dimension (solution vector lenght)
		long double linf, lsup;											// lower and upper limits for the real numbers
		int change_type;												// change type: 0: no change 1:type 1.1 , 2:type 1.2 , 3:type 1.3 , 4:type 2.1 , 5:type 2.2 , 6. type 2.3, 7:type 3.1
		int os;															// order of the template (schema)
		int rho0_flag;													// flag indicating that rho=0
		int *b_t12, *l_t21,  *s; 										// control parameters
		long double *m_t22, *m_t23, **M_t11, **M_t13;					// control parameters
		long double *xb, df_t31, dx_range;								// control parameters
		void scaling( long double *x , long double *xs );
		void unscaling( long double *x  );
	public:		
		dop(int l_p, long double linf_p,  long double lsup_p);		
	    ~dop(void);
		long double transform( long double *x ,  long double *x_new);						// fitness transformation
		void change( int change_type_p, long double rho, long double f_range, long double *x_template ); 	// change in the DOP
}; 



/******************************************************************************\
*								Constructor													   *
\******************************************************************************/
dop::dop(int l_p, long double linf_p,  long double lsup_p){
	
	// Parameters
	l=l_p;				
	linf=linf_p;
	lsup=lsup_p;
	//	dx_range=0.01*fabsl(lsup-linf)/sqrt(l);		// parameter used to define the range close to the schema variable
	dx_range=0.05*fabsl(lsup-linf);		// parameter used to define the range close to the schema variable

	// Memory allocation 
	M_t11=aloc_matrixd(l,l);
	b_t12=aloc_vectori(l);
	M_t13=aloc_matrixd(l,l);
	l_t21=aloc_vectori(l);
	m_t22=aloc_vectord(l);
	m_t23=aloc_vectord(l);
	s=aloc_vectori(l);
	xb=aloc_vectord(l);		// xbest in previous environment
	
	change_type=0;		// stationary landscape
		
}


/******************************************************************************\
*								 Destructor													   *
\******************************************************************************/
dop::~dop(void){

	desaloc_matrixd(M_t11,l);
	delete [] b_t12;
	desaloc_matrixd(M_t13,l);
	delete [] l_t21;
	delete [] m_t22;
	delete [] m_t23;
	delete [] s;
	delete [] xb;

	               
}  


/******************************************************************************\
*			Scaling Between -1 and 1 (used by DOP type 1.1)		 			   *
\******************************************************************************/
void dop::scaling( long double *x , long double *xs ){
	int i;
	
	for (i=0;i<l;i++)
		xs[i] = (long double) 2.0*(x[i]-linf)/(lsup-linf) - 1.0;
		
} 


/******************************************************************************\
*			Unscaling  (used by DOP type 1.1)		 			   			   *
\******************************************************************************/
void dop::unscaling( long double *x  ){
	int i;
	
	for (i=0;i<l;i++){
		x[i] = (long double) (x[i]+1.0)*(lsup-linf)/2.0 + linf;
		if ( x[i]> lsup ) 
			x[i] = lsup;
		else if ( x[i]< linf)  
			x[i] = linf;					
	}
	
} 


/******************************************************************************\
*								 Fitness transformation			 			   *
\******************************************************************************/
long double dop::transform( long double *x ,  long double *x_new){
	int i, schema_flag;
	long double df, *xs;
	
	if (change_type==0){
		// No change
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;		
	} 
	else if (change_type==1){
		// DOP Type 1.1 (DOPs with rotation (see Ref. 2)	)
		xs=aloc_vectord(l);
		scaling(x, xs);	
		multMatrixVect(x_new, M_t11, l, l, xs, l); 				// Eq. 4 of ref. 2
		delete [] xs;
		unscaling(x_new);
		df=0.0;
	} 
	else if (change_type==2){
		// DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)
		// Eq. 37: observation - we don´t use a matrix B_t12 because of it costs O(l^2)
		for (i=0;i<l;i++)
			x_new[i]=x[b_t12[i]];
		df=0.0;
	} 
	else if (change_type==3){
		// DOP Type 1.3 (DOPs with rotation of the candidate solutions according to a template)	
		// Comparing x with the template and changing the respective variables
        schema_flag=1;           // 1 if x in s and 0 otherwise
        if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
            if (s[i]>=0)
                if ( fabsl(x[i]-xb[i])>dx_range )
                    schema_flag=0;
            i++;
        }
        // Eq. 39
        if (schema_flag==1){
        	xs=aloc_vectord(l);
			scaling(x, xs);	
			multMatrixVect(x_new, M_t13, l, l, xs, l); 				// Eq. 4 of ref. 2
			delete [] xs;
			unscaling(x_new);
			df=0.0;	        	
		}
		else {		
			for (i=0;i<l;i++)
				x_new[i]=x[i];	
		}                      			
		df=0.0;
	}
	else if (change_type==4){
		// DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)	
		// Eq. 43: observation - we don´t use a matrix L_t21 because of it costs O(l^2)
		for (i=0;i<l;i++)
			x_new[i]=x[l_t21[i]];
		df=0.0;
	} 
	else if (change_type==5){
	   // DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates: m is random)															
	   // Comparing x with the template and changing the respective variables
        schema_flag=1;           // 1 if x in s and 0 otherwise
        if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
            if (s[i]>=0)
                if ( fabsl(x[i]-xb[i])>dx_range )
                    schema_flag=0;
            i++;
        }
        // Eq. 45
        if (schema_flag==1){        
			for (i=0;i<l;i++)
				x_new[i]=m_t22[i];
		}
		else {		
			for (i=0;i<l;i++)
				x_new[i]=x[i];	
		}			
		df=0.0;
	} 
	else if (change_type==6){
		// DOP Type 2.3 (DOPs obtained by copying decision variables according to a set of templates: m=xb)															
		// Comparing x with the template and changing the respective variables
        schema_flag=1;           // 1 if x in s and 0 otherwise
        if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
            if (s[i]>=0)
                if ( fabsl(x[i]-xb[i])>dx_range )
                    schema_flag=0;
            i++;
        }
        // Eq. 45
        if (schema_flag==1){        
			for (i=0;i<l;i++)
				x_new[i]=m_t23[i];
		}
		else {		
			for (i=0;i<l;i++)
				x_new[i]=x[i];	
		}			
		df=0.0;
	} 
	else if (change_type==7){
		// DOP Type 3.1: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;
		// Comparing x with the template and changing the respective variables
        schema_flag=1;           // 1 if x in s and 0 otherwise
        if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
            if (s[i]>=0)
                if ( fabsl(x[i]-xb[i])>dx_range )
                    schema_flag=0;
            i++;
        }
        // Eq. 47, 48, 49
        if (schema_flag==1)
            df=df+df_t31;
                  	
	} 
	
	return df;
}


/******************************************************************************\
*								Change in the DOP							   *
\******************************************************************************/
void dop::change( int change_type_p, long double rho, long double f_range, long double *x_template ){

	int i, j, temp, line_1, line_2, *perm_v, *temp_v;
	long double df_range, theta, cos_theta, sin_theta;
	
	change_type=change_type_p;

	// random permutation of a vector of integers
	perm_v=aloc_vectori(l);			
	temp_v=aloc_vectori(l);			
	for (i=0;i<l;i++)
		temp_v[i]=i;	
 	
 	// order of the schema (used in types 1.3, 2.2., 2.3 and 3.1)
 	if (rho>0.0){
    	os=floorl( -log2l(rho/2.0) );  
    	if (os>l)
    		os=l; 
    	//os=ceill( (1.0-rho)*l );   // obs.: depends on l
    	rho0_flag=0;
	}
	else{
		os=l;  
		rho0_flag=1;
	}
			
	if (change_type==1){
		// DOP Type 1.1 (DOPs with rotation)	see ref. 2
		// See Alg. 1 - ref. 2
		// Identity Matrix
		for (i=0;i<l;i++) {
			M_t11[i][i]=1.0;
			for (j=i+1;j<l;j++){
				M_t11[i][j]=0.0;	
				M_t11[j][i]=0.0;
			}
		}	
		// Angle
		theta = rho*PI;
		cos_theta = (long double) cos(theta);
		sin_theta = (long double) sin(theta);
		// rotations
		rand_perm(temp_v,perm_v,l);
		for (i=0;i<l-1;i=i+2) {
			M_t11[ perm_v[i] ][ perm_v[i] ]     = cos_theta;
			M_t11[ perm_v[i] ][ perm_v[i+1] ]   = -sin_theta;
			M_t11[ perm_v[i+1] ][ perm_v[i] ]   = sin_theta;
			M_t11[ perm_v[i+1] ][ perm_v[i+1] ] = cos_theta;
		}		
	} 
	else if (change_type==2){
		// DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation matrix)
		// permuting 2 lines of matrix B_t12
		// Eq. 38
		for (i=0;i<l;i++)
			b_t12[i]=i;
		for (i=0;i<rho*l;i++){		
			line_1=random_int(0,l-1);
			line_2=random_int(0,l-1);
			while (line_2==line_1)
				line_2=random_int(0,l-1);							
			temp=b_t12[line_1];			
			b_t12[line_1]=b_t12[line_2];
			b_t12[line_2]=temp;			
		}
	} 
	else if (change_type==3){
		// DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a set of templates)
		// Eq. 40	
		rand_perm(temp_v,perm_v,l);
		i=0;
		while (i<l){
			xb[i]=x_template[i];		
			if(i<os)	
				s[perm_v[i]]=1;					// fixed variables
			else 					
				s[perm_v[i]]=-1;				// -1 indicates don´t variables
			i++;
		}
		// See Alg. 1 - ref. 2
		// Identity Matrix
		for (i=0;i<l;i++) {
			M_t13[i][i]=1.0;
			for (j=i+1;j<l;j++){
				M_t13[i][j]=0.0;	
				M_t13[j][i]=0.0;
			}
		}	
		// Angle
		theta = rho*PI;
		cos_theta = (long double) cos(theta);
		sin_theta = (long double) sin(theta);
		// rotations
		rand_perm(temp_v,perm_v,l);
		for (i=0;i<l-1;i=i+2) {
			M_t13[ perm_v[i] ][ perm_v[i] ]     = cos_theta;
			M_t13[ perm_v[i] ][ perm_v[i+1] ]   = -sin_theta;
			M_t13[ perm_v[i+1] ][ perm_v[i] ]   = sin_theta;
			M_t13[ perm_v[i+1] ][ perm_v[i+1] ] = cos_theta;
		}					
	}
	else if (change_type==4){
		// DOP Type 2.1 (DOPs obtained by copying elements of the decision variables according to a linear transformation)	
		// Eq. 44
		for (i=0;i<l;i++)
			l_t21[i]=i;
		temp=(int) ceil(rho*l/2.0);
		for (i=0;i<temp;i++){		
			line_1=random_int(0,l-1);
			line_2=random_int(0,l-1);
			while (line_2==line_1)
				line_2=random_int(0,l-1);	
			l_t21[line_1]=l_t21[line_2];
		}
	} 
	else if (change_type==5){
		// DOP Type 2.2 (DOPs obtained by copying decision variables according to a set of templates: m random)
		// Eq. 46
		rand_perm(temp_v,perm_v,l);
		i=0;
		while (i<l){
			xb[i]=x_template[i];		
			if(i<os){			
				s[perm_v[i]]=1;					// fixed variables
				m_t22[perm_v[i]]=x_template[perm_v[i]];
			}
			else{						
				s[perm_v[i]]=-1;				// -1 indicates don´t care variables
				m_t22[perm_v[i]]=(long double) ( random_dou() )*(lsup-linf)+linf;	 // random between lim_inf and lim_sup
			}
			i++;
		}
	
	} 
	else if (change_type==6){
		// DOP Type 2.3 (DOPs obtained by copying decision variables according to a set of templates: m=xb)
		// Eq. 46															
		rand_perm(temp_v,perm_v,l);
		i=0;
		while (i<l){
			xb[i]=x_template[i];		
			if(i<os){			
				s[perm_v[i]]=1;					// fixed variables
				m_t23[perm_v[i]]=x_template[perm_v[i]];
			}
			else{						
				s[perm_v[i]]=-1;				// -1 indicates don´t care variables
				m_t23[perm_v[i]]=x_template[perm_v[i]];	 // random between lim_inf and lim_sup
			}
			i++;
		}				
	} 
	else if (change_type==7){
		// DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	
		// Eq. 49			
		df_range=f_range*rho;		
		rand_perm(temp_v,perm_v,l);
		i=0;
		while (i<l){		
			if(i<os)			
				s[perm_v[i]]=1;
			else 					
				s[perm_v[i]]=-1;				// -1 indicates don´t care bits					
			i++;
		}
		df_t31=random_dou()*(2.0*df_range)-df_range;			// observation: here, uniform distribution is used instead of normal distribution
		
	}
	
	
	delete [] temp_v;
	delete [] perm_v;
}


