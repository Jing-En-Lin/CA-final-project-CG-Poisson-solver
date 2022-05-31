#include <cstdio>
#include <cstdlib>
#include <omp.h>
#include <iostream>
#include <cmath>

using namespace std;

//constants
   double L   = 1.0;   // 1-D computational domain size
   const int k   = 256;   // number of equally spaced sampling points on both x and y direction
   const int N = k*k;
   const int g = k-2;
   const int NThread = 5;
   double epsilon = 5.0e-7; //acceptable final error 
   double amp = 1.0;
   int totalStepslimit = 2000;

   // initial condition
   int    steps = 0;                     // number of iterations
   double phi[N]={0};
   double rho[N];
   double phi_ref[N]={0};
   double err = 0;
   double dx;
   double dA ;  
   double phi_inner[g*g]={0};
   double rho_inner[g*g];
   double r[g*g];
   double d[g*g];
   double delta_new = 0;
   double delta_final = 0;
   double delta_old = 0;
   double q[g*g]={0};
   double alpha;
   double beta;
   int tid;
   


double rho_ref_func(double amp,double x, double y ){
	   return amp*(( y + 1.0)*( 4.0*x*x*x + 4.0*x*x - 6.0*x - 2.0 ) + ( x + 1.0)*( 4.0*y*y*y + 4.0*y*y - 6.0*y - 2.0 ))*exp(-1.0*y*y-1.0*x*x);	
	}

   
// can be valid only on boundary
double phi_ref_func( double amp, double x, double y ){
       return amp*( x + 1 )*( y + 1 )*exp( -1*x*x-y*y );	
    }
    // setting corners twice, which is useless but easy to understand
void   phi_BC(double phi[N]){
	    //top boundary 
	    //#pragma omp parallel
	    {
		  //#pragma omp for
	      for(int i=0;i<k;i++)
	       {
		    phi[i]=phi_ref_func( amp, i*dx, L);
		    phi[i*k]=phi_ref_func( amp, 0, L-i*dx);
		    phi[i*k+k-1]=phi_ref_func( amp, L , L-i*dx);
		    phi[i+(k-1)*k]=phi_ref_func( amp, (i%k)*dx, 0);
		   }
		 
	   }
		  
//	   // i = j*k  left
//	   for(int j=0;j<k;j++)
//	     phi[j*k]=phi_ref_func( amp, 0, L-j*dx);
//	   // i = j*k+k-1  right
//	   for(int j=0;j<k;j++)
//	     phi[j*k+k-1]=phi_ref_func( amp, L , L-j*dx);
//	     // botton
//	   for(int i=(k-1)*k;i<N;i++)
//	     phi[i]=phi_ref_func( amp, (i%k)*dx, 0);
    }    
   
void   rho_BC(double rho[N]){
	    // setting corner only once since we use the same array 
	   rho[k+1] =  ( phi_ref[k]+phi_ref[1] )/dA + rho[k+1];
	   rho[2*k-2]= ( phi_ref[k-2]+phi_ref[2*k-1] )/dA + rho[2*k-2];
	   rho[(k-2)*k+1] = (phi_ref[(k-1)*k+1]+phi_ref[(k-2)*k] )/dA + rho[(k-2)*k+1];
	   rho[(k-1)*k-2] = (phi_ref[(k-1)*k-1]+phi_ref[k*k-2] )/dA + rho[(k-1)*k-2];
	   # pragma omp parallel
	   {
	     #pragma omp for
	     for(int i=0;i<k-4;i++)
	     {
	     	rho[i+k+2]= phi_ref[i+2]/dA + rho[i+k+2];
	    	rho[(i+2)*k+1]= phi_ref[(i+2)*k]/dA + rho[(i+2)*k+1];
	    	rho[(i+3)*k-2]= phi_ref[(i+3)*k-1]/dA + rho[(i+3)*k-2];
	    	rho[i+(k-2)*k+2]= phi_ref[i+(k-1)*k+2]/dA + rho[i+(k-2)*k+2];
	     }
       }
//	   for(int i=k+2;i<2*k-2;i++)
//	   	rho[i]= phi_ref[i-k]/dA + rho[i];
//	    
//	   for(int j=2;j<k-2;j++)
//	   { int i=j*k+1;
//	   	rho[i]= phi_ref[i-1]/dA + rho[i];
//	   }
//	   
//	   for(int j=3;j<k-1;j++)
//	   {
//	    int i = j*k-2;
//	   	rho[i]= phi_ref[i+1]/dA + rho[i];
//	   }
//	   
//	   for(int i=(k-2)*k+2;i<(k-1)*k-2;i++)
//	    rho[i]= phi_ref[i+k]/dA + rho[i];
	   	
	   //barrier
	   
    } 
    // it will give us the inner part of x 
void   extract(double x_ext[g*g], double x[N] ){
	   for(int i=1;i<k-1;i++)
	    for(int j=1;j<k-1;j++)
	      	x_ext[i-1+(j-1)*g]=x[i+j*k];
	}

void   merge(double x_ext[g*g], double x[N]){
	   
	   for(int i=1;i<k-1;i++)
	    for(int j=1;j<k-1;j++)
	      	x[i+j*k] = x_ext[i-1+(j-1)*g];
	
    }
	
	  
//  r=b-Ax   matrix A belongs to M_g^2 * g^2       
void   matrix_mult(double r[g*g],double b[g*g],double x[g*g]){
	   //calculate corner first
	   r[0] = b[0] + ( -4*x[0] +x[1] +x[g] )/dA;
	   r[g-1] = b[g-1] + ( -4*x[g-1] +x[g-2] +x[2*g-1] )/dA;
	   r[(g-1)*g] = b[(g-1)*g] + ( -4*x[(g-1)*g] +x[(g-1)*g+1] +x[(g-2)*g] )/dA;
	   r[g*g-1] = b[g*g-1] + ( -4*x[g*g-1] +x[g*g-2] +x[(g-1)*g-1] )/dA;
	   
	   # pragma omp parallel
    {
	   //calculate side next
	   # pragma omp for nowait
	   for(int i=1;i<g-1;i++)
	      {
		  r[i] = b[i] + ( -4*x[i] +x[i-1] +x[i+1] +x[i+g] )/dA;
	      r[i+(g-1)*g] = b[i+(g-1)*g] + ( -4*x[i+(g-1)*g] +x[i+(g-1)*g-1] +x[i+(g-1)*g+1] +x[i+(g-1)*g-g] )/dA;
	      r[i*g] = b[i*g] + ( -4*x[i*g] +x[i*g+1] +x[(i-1)*g] +x[(i+1)*g] )/dA;
	      r[(i+1)*g-1] = b[(i+1)*g-1] + ( -4*x[(i+1)*g-1] +x[(i+1)*g-2] +x[i*g-1] +x[(i+2)*g-1] )/dA;
         }
	   
	   // i+(g-1)*g	
//	   # pragma omp for nowait
//	   for(int i=(g-1)*g+1;i<g*g-1;i++)
//	      r[i] = b[i] + ( -4*x[i] +x[i-1] +x[i+1] +x[i-g] )/dA;
	   	
	   
	   // i = j*g
//	   # pragma omp for nowait
//	   for(int j=1;j<g-1;j++)
//	      r[j*g] = b[j*g] + ( -4*x[j*g] +x[j*g+1] +x[(j-1)*g] +x[(j+1)*g] )/dA;
	   	
	   
	   // i = j*g-1     //(i+1)
//	   # pragma omp for nowait
//	   for(int j=2;j<g;j++)
//	      r[j*g-1] = b[j*g-1] + ( -4*x[j*g-1] +x[j*g-2] +x[(j-1)*g-1] +x[(j+1)*g-1] )/dA;
	   	
	   
	   //inner part   //remember to collapse
	   # pragma omp for  collapse(2)
	   for(int i=1;i<g-1;i++)
	   {
	   	 for(int j=1;j<g-1;j++)
	   	 {
	   	 	r[i+j*g]= b[i+j*g] + ( - 4*x[i+j*g] + x[i+j*g-1] + x[i+j*g+1] + x[i+(j-1)*g] + x[i+(j+1)*g] )/dA;
		 }    
	   }
    }
}    

//  q=Ad    
void   matrix_mult2(double q[g*g],double d[g*g]){
	   //calculate corner first
	   q[0] = ( 4*d[0] -d[1] -d[g] )/dA;
	   q[g-1] = ( 4*d[g-1] -d[g-2] -d[2*g-1] )/dA;
	   q[(g-1)*g] = ( 4*d[(g-1)*g] -d[(g-1)*g+1] -d[(g-2)*g] )/dA;
	   q[g*g-1] = ( 4*d[g*g-1] -d[g*g-2] -d[(g-1)*g-1] )/dA;
	   
	   # pragma omp parallel
    {
	   //calculate side next
	   # pragma omp for nowait
	   for(int i=1;i<g-1;i++)
	      {
		  q[i] = ( 4*d[i] -d[i-1] -d[i+1] -d[i+g] )/dA;
		  q[i+(g-1)*g] = ( 4*d[i+(g-1)*g] -d[i+(g-1)*g-1] -d[i+(g-1)*g+1] -d[i+(g-1)*g-g] )/dA;
		  q[i*g] = ( 4*d[i*g] -d[i*g+1] -d[(i-1)*g] -d[(i+1)*g] )/dA;
		  q[(i+1)*g-1] = ( 4*d[(i+1)*g-1] -d[(i+1)*g-2] -d[i*g-1] -d[(i+2)*g-1] )/dA;
	   	  } 
	   	  
	   	  
//	   # pragma omp for nowait
//	   for(int i=(g-1)*g+1;i<g*g-1;i++)
//	      q[i] = ( 4*d[i] -d[i-1] -d[i+1] -d[i-g] )/dA;
	   	
	   
	   // i = j*g
//	   # pragma omp for nowait
//	   for(int j=1;j<g-1;j++)
//	      q[j*g] = ( 4*d[j*g] -d[j*g+1] -d[(j-1)*g] -d[(j+1)*g] )/dA;
	   	
	   
	   // i = j*g-1   
//	   # pragma omp for nowait
//	   for(int j=2;j<g;j++)
//	      {q[j*g-1] = ( 4*d[j*g-1] -d[j*g-2] -d[(j-1)*g-1] -d[(j+1)*g-1] )/dA;
//	   	     if(tid == 5) 
//			 printf("%d",tid);
//			 } 
	   
	   //inner part   //remember to collapse
	   # pragma omp for collapse(2)
	   for(int i=1;i<g-1;i++)
	   {
	   	 for(int j=1;j<g-1;j++)
	   	 {
	   	 	q[i+j*g]= ( 4*d[i+j*g] - d[i+j*g-1] - d[i+j*g+1] - d[i+(j-1)*g] - d[i+(j+1)*g] )/dA;
	   	 	
		 }    
	   }
    }
}    


//innner product = d^T*q
double inn_product(double d[g*g],double q[g*g]){
	   double temp = 0;
	   for(int i=0;i<g*g;i++)
	   {
	   	temp += q[i]*d[i];
	   }
	 return temp;
    }






 int main( int argc, char *argv[] )
 {
 	
 	//set the thread number 
 	omp_set_num_threads( 6 );
 	
 	
   // derived constants 
  dx      = L/(k-1);    //spatial resolution
  dA      = dx*dx;
     
    
	
	//initialize rho to reference function 
	 for(int i=0;i<N;i++)
    { //rho must have negative sign so that A can be positive-definite
	  rho[i] = -rho_ref_func( amp, (i%k)*dx, L-(i/k)*dx );
	  phi_ref[i] = phi_ref_func( amp, (i%k)*dx, L-(i/k)*dx );
 	}
// 	for(int i=0;i<N;i++)
// 	{
// 		printf("%10.3e  ",rho[i]);
// 		if(i%k==k-1)
// 		printf("\n");
//	 }
// 	printf("\n\n\n");
	
 	
    phi_BC(phi);   
	// use rho to force phi maintaining boundary condition 
	rho_BC(rho);
	
//	for(int i=0;i<N;i++)
// 	{
// 		printf("%10.3e  ",rho[i]);
// 		if(i%k==k-1)
// 		printf("\n");
//	 }
	
	// extract the inner part which we need solve
	 
	extract(rho_inner,rho); 
    
    

     
// The CG iteration part (start to parellel)
    
	 //r=b-Ax 
	 // In fact, phi_inner is zero so r = b
	 matrix_mult(r,rho_inner,phi_inner);
	 //d=r
     for(int i=0;i<g*g;i++)
     {
        d[i] = r[i];	
     }
     //delta_new = <r,r>
     delta_new = inn_product(r,r);
    
     delta_final = epsilon*epsilon*delta_new;
    
	//loop
    for(int stepNumber= 0; stepNumber< totalStepslimit; stepNumber++ )     // stepNumber is "steps"
   {
     
       
       
       //tid = omp_get_thread_num();	
       
       
	   matrix_mult2(q,d);
     
       alpha = delta_new/inn_product(d,q);
       
	   #pragma omp parallel
	   {
		  
       #pragma omp for
       for(int i=0;i<g*g;i++)
         phi_inner[i] = phi_inner[i] + alpha*d[i];
        
       }
        
	   // this correction is used only in the very acurate demand
	   if( stepNumber == g )
          matrix_mult(r,rho_inner,phi_inner);

	   else
	     {
	       #pragma omp parallel
	       {
		   #pragma omp for
	       for(int i=0;i<g*g;i++)
	          r[i] = r[i]-alpha*q[i];
           }
         }
	
	   delta_old = delta_new;
	   delta_new = inn_product(r,r);
	   beta = delta_new/delta_old;
	
	   #pragma omp parallel
	 {
	 	#pragma omp for
	    for(int i=0;i<g*g;i++)
          d[i]=r[i]+beta*d[i];  
		  
     }  	
	 
       merge(phi_inner,phi); 	
	    
     
    err = 0;
    for(int i=0;i<N;i++)
    err+= abs((phi[i]-phi_ref[i])/phi_ref[i]);
    
    err/=float (g*g);          //normalization

    steps++;       //  update steps
    if(steps %100 == 0) 
	printf("steps = %d  err= %10.3e \n",steps, err);
    
	if(delta_new<delta_final)
	   break;

   }
   
 	
  	

//	for(int i=0;i<N;i++)
// 	{
// 		printf("%10.3e  ",phi[i]);
// 		if(i%k==k-1)
// 		printf("\n");
//    }

  printf("steps = %d  fin_err= %10.3e \n",steps, err);
  return EXIT_SUCCESS;
 }
	
	
