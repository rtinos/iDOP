/*****************************************************************************************************************************************************\
*								 Class: DOP																											  *	
* Main reference:																																	  *
* Tinos, R. & Yang, S. (2019). "A Framework for Inducing Artificial Changes in Optimization Problems", Submitted to Information Sciences.			  *
* Aditional references:			 				   																									  *
*  Ref. 1: Tinos, R. & Yang, S. (2014). "Analysis of fitness landscape modifications in evolutionary dynamic optimization", Information Sciences, 282.*
*  Ref. 2: TinosR. & Yang, S. (2007). "Continuous dynamic problem generators for evolutionary algorithms", Proc. of IEEE CEC, 236-243.				  *
* 																																					  *
* DOP types (Pseudo-Boolean):																														  *
* 	DOP Type 1: DOP with permutation 																												  *
* 		DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)						  *
* 		DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation of x)		  *
* 		DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to template)				  *
* 	DOP Type 2: Single time-dependent DOP obtained by copy of decision variables																	  *
* 		DOP Type 2.1 (DOPs obtained by copying elements of the decision variables in vector x)														  *
* 		DOP Type 2.2 (DOPs obtained by copying decision variables from a random solution according to a template)									  *
* 		DOP Type 2.3 (DOPs obtained by copying decision variables from the current best solution according to a template)							  *
* 	DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a template													  *	
* 																																					  *
* DOP types (Real):																																	  *
* 	DOP Type 1: DOP with rotation 																												      *
* 		DOP Type 1.1 (DOPs with rotation of the candidate solutions defined by candidate solution exchanges of the XOR type)						  *
* 		DOP Type 1.2 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a permutation of x)		  *
* 		DOP Type 1.3 (DOPs with rotation of the candidate solutions defined by decision variable exchanges according to template)					  *
* 	DOP Type 2: Single time-dependent DOP obtained by copy of decision variables																	  *
* 		DOP Type 2.1 (DOPs obtained by copying elements of the decision variables in vector x)														  *
* 		DOP Type 2.2 (DOPs obtained by copying decision variables from a random solution according to a template)									  *
* 		DOP Type 2.3 (DOPs obtained by copying decision variables from the current best solution according to a template)							  *
* 	DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a template													  *		
\*****************************************************************************************************************************************************/


#include <cmath>
#include <cstdlib>
#define N_max_s 3

class dop{
	private:	
		int l;															// size of the chromosome
		int change_type;												// change type: 0: no change 1:type 1.1 , 2:type 1.2 , 3:type 1.3 , 4:type 2.1 , 5:type 2.2 , 6. type 2.3, 7:type 3.1
		int os;															// order of the templates (schemata)
		int rho0_flag;													// flag indicating that rho=0
		int *m_t11, *b_t12, *m_t13, *l_t21, *m_t22, *m_t23, *s; 		// control parameters
		double df_t30;													// control parameters
	public:		
		dop(int l_p);		
	    ~dop(void);
		double transform( int *x ,  int *x_new);						// fitness transformation
		void change( int change_type_p, double rho, double f_range, int *x_template ); 	// change in the DOP
}; 



/******************************************************************************\
*								Constructor													   *
\******************************************************************************/
dop::dop(int l_p){
	
	// Parameters
	l=l_p;				// dimension (solution vector lenght)
	
	// Memory allocation 
	m_t11=aloc_vectori(l);
	b_t12=aloc_vectori(l);
	m_t13=aloc_vectori(l);
	l_t21=aloc_vectori(l);
	m_t22=aloc_vectori(l);
	m_t23=aloc_vectori(l);
	s=aloc_vectori(l);
	
	change_type=0;		// stationary landscape
		
}


/******************************************************************************\
*								 Destructor													   *
\******************************************************************************/
dop::~dop(void){

	delete [] m_t11;
	delete [] b_t12;
	delete [] m_t13;
	delete [] l_t21;
	delete [] m_t22;
	delete [] m_t23;
	delete [] s;

	               
}  


/******************************************************************************\
*								 Fitness transformation			 			   *
\******************************************************************************/
double dop::transform( int *x ,  int *x_new){
	int i, schema_flag;
	double df;
	
	if (change_type==0){
		// No change
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;		
	} 
	else if (change_type==1){
		// DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)	
		XOR(x,m_t11,x_new,l);						// Eq. 35
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
		// DOP Type 1.3 (DOPs with permutation of the candidate solutions defined by decision variable exchanges according to a template)	
		// Comparing x with the template and changing the respective bits
        schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
         if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
            if (s[i]>=0)
                if (s[i] != x[i] )
                    schema_flag=0;
            i++;
        }
        // Eq. 39
        if (schema_flag==1)
			XOR(x,m_t13,x_new,l);  
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
		// Comparing x with the template and changing the respective bits
	    schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
	    if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
        	if (s[i]>=0)
            	if (s[i] != x[i] )
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
		// Comparing x with the template and changing the respective bits
	     schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
	    if (rho0_flag==1)
        	schema_flag=0;
        i=0;
        while (schema_flag>0 && i<l){
        	if (s[i]>=0)
            	if (s[i] != x[i] )
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
		// DOP Type 3: Single time-dependent DOPs obtained by adding fitness terms according to a set of templates	
		for (i=0;i<l;i++)
			x_new[i]=x[i];
		df=0.0;
		// Comparing x with the  template and changing the respective bits
        schema_flag=1 ;           // 0 if schema m is not present in x and 1 otherwise
        if (rho0_flag==1)
        	schema_flag=0;
	    i=0;
	    while (schema_flag>0 && i<l){
	        if (s[i]>=0)
	            if (s[i] != x[i] )
	               	schema_flag=0;
	        i++;
	    }
        // Eq. 47, 48, 49
        if (schema_flag==1)
            df=df+df_t30;
                  	
	} 
	
	return df;
}


/******************************************************************************\
*								Change in the DOP							   *
\******************************************************************************/
void dop::change( int change_type_p, double rho, double f_range, int *x_template ){

	int i, l_minus_order2, temp, line_1, line_2, *perm_v, *temp_v, *r;
	double df_range;
	
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
   
    	rho0_flag=0;
	}
	else{
		os=l;  
		rho0_flag=1;
	}
	 l_minus_order2 = (l-os)/2;          // the order of a schema is the number of fixed positions in the schema 	
	 	
	if (change_type==1){
		// DOP Type 1.1 (DOPs with permutation of the candidate solutions defined by candidate solution exchanges of the XOR type)	
		rand_perm(temp_v,perm_v,l);
		r=aloc_vectori(l);
		i=0;
		while (i<l){	
			m_t11[i]=0;		// here, m_t11 and r will be equal; we are keeping both variable to keep the same structure of the DOP generator 
			if(i<rho*l)
				r[perm_v[i]]=1;
			else
				r[perm_v[i]]=0;
			i++;
		}
		for (i=0;i<l;i++)
			temp_v[i]=m_t11[i];
		XOR(temp_v,r,m_t11,l);		// Eq. 36		
		delete [] r;	
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
			if(i<os){	
				s[perm_v[i]]=x_template[perm_v[i]];				
				m_t13[perm_v[i]]=0;
			}
			else {					
				s[perm_v[i]]=-1;				// -1 indicates don´t care bits
				if(i>l_minus_order2)
					m_t13[perm_v[i]]=1;
				else
					m_t13[perm_v[i]]=0;
			}
			i++;
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
			if(i<os){				
				s[perm_v[i]]=x_template[perm_v[i]];
				m_t22[perm_v[i]]=s[perm_v[i]];
			}
			else {					
				s[perm_v[i]]=-1;				// -1 indicates don´t care bits
				m_t22[perm_v[i]]=random_int(0,1);
					
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
			if(i<os){				
				s[perm_v[i]]=x_template[perm_v[i]];
				m_t23[perm_v[i]]=x_template[perm_v[i]];
			}
			else {					
				s[perm_v[i]]=-1;				// -1 indicates don´t care bits
				m_t23[perm_v[i]]=x_template[perm_v[i]];
					
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
				s[perm_v[i]]=x_template[perm_v[i]];
			else 					
				s[perm_v[i]]=-1;				// -1 indicates don´t care bits					
			i++;
		}
		df_t30=random_dou()*(2.0*df_range)-df_range;			// observation: here, uniform distribution is used instead of normal distribution
		
	}
	
	
	delete [] temp_v;
	delete [] perm_v;
}


