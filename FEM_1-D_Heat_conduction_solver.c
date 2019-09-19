
#include <math.h>
#include <stdio.h>

int main (void)
{
	
	int   m = 2,n = 2, p, c,d, i,j,sum,o;
	float kxAl,kxCu , length ,dia, kalpre, kcupre , area, deltax,e ;
	float kal[10][10], kcu[10][10] , kas1[10][10], kas2[10][10],kas[10][10] ;
	float Diaginv[10][10], E[10][10],  F [10][10] , RHS1 [10][1] , RHS2[10][1] ;
	float k[2][2];
    float q[10][1];
    float approxT[10][1];
	
	
	printf("********** 1-D heat Conduction ***** FEM Solver ********** \n");
	
	printf("\n");
	
	printf("Problem Statement : \n");
	
	printf(" Considering a Hot rod of composite material like Al and Cu of equal length ");
	printf("\n");
	
	printf(" Input - Heat Flux and properties");
	printf("\n");
	
	printf(" Output - we will obtain the temperature distribution");
	printf("\n");
	
	printf("Enter the diameter of the rod in meters \n");
	scanf("%f", &dia);
	printf("\n");
	area =  (3.14*dia*dia)/4 ; 
    printf("Area of the rod: %f \n", area);
    printf("\n");
	printf("Only Even number elements division to be entered \n");
	printf("\n");
	printf("Enter the no of elements need to divided (Consider 4 elements for Convenience) \n");
	scanf("%f", &e);
	printf("\n");
	printf("Enter the length of the rod in meters \n");
	scanf("%f", &length);
	printf("\n");
	deltax = (length)/e ;
	
	printf("Enter the Thermal conductivity of Aluminium \n");
	scanf("%f", & kxAl);
	
	printf("\n");
	printf("Enter the Thermal conductivity of copper \n");
	scanf("%f", & kxCu);
	printf("\n");
//******************************  LHS     ********************************
	
    printf("Grid size (Delta X) : %f \n", deltax);
	kalpre = (kxAl*area)/deltax ;
//	printf("Kal pre  : %f \n", kalpre);
	// Stiffness matrix for aluminium
	for(p=0;p<(e/2);p++)
   {
   printf("Enter the elements of local conductance matrices for aluminium \n");
 
   for (c = p; c < m+p; c++){
   	
   	  for (d = p; d < n+p; d++){
   	  	
         scanf("%f", &k[c][d]);
		 }
   }
   

    for (c=p ; c<m+p; c++){
 	for (d =p; d<n+p;d++){
 		
 		
        kal[c][d] = kalpre*k[c][d];
 		printf("%f\t", kal[c][d]);
 }
  printf("\n");
}
}
	// assembling of conductance matricesof Aluminium
  printf("The global conductance matrices for aluminum : \n");

    for (c=0 ; c<(e/2)+1; c++){
 	for (d=0;  d<(e/2)+1; d++){
 		
 		if (c==d && c!=0 && d!=0 && c!=e/2 && d!= e/2) {
 			
 			kas1[c][d] = kal[c][d] + kal[c][d];	
		 }
 		else	
 		kas1[c][d] =  kal[c][d]  ; 	
		
 		printf("%f\t", kas1[c][d]);
 		
 }
 	 
	  printf("\n");
}
	// conductance matrices for copper
	kcupre = (kxCu*area)/deltax ;
	
	for(p=(e/2);p<e ;p++)
   {
 printf("Enter the elements of local conductance matrices  for copper\n");
 
   for (c = p; c < m+p; c++)
      for (d = p; d < n+p; d++)
         scanf("%f", &k[c][d]);
 
 
    for (c=p ; c<m+p; c++){
 	for (d=p; d<n+p;d++){
 		
 		
        kcu[c][d] = kcupre *k[c][d];
 		printf("%f\t", kcu[c][d]);
 }
  printf("\n");
}
}

// assembling of conductance matrices of copper
 

  printf("The global conductance matrices for copper : \n");

    for (c=(e/2); c<e+1; c++){
 	for (d=(e/2); d<e+1; d++){
 		
 		if (c==d && c!= (e/2) && d!=(e/2) && c!=e && d!= e) {
 			
 			kas2[c][d] = kcu[c][d] + kcu[c][d];	
		 }
 		else	
 		kas2[c][d] =  kcu[c][d]  ; 	
		
 		printf("%f\t", kas2[c][d]);
 		
 }
 	 
	  printf("\n");
}

// complete assembling
 printf("The global conductance matrices : \n");

    for (c=0; c<e+1; c++){
 	for (d=0; d<e+1; d++){
 		
 	//	if (c==d && c!=0 && d!=0 && c!=e && d!= e) {
 		if (c==d)	{	 
 			kas[c][d] = kas1[c][d] + kas2[c][d];	
		 }
 		else	
 		kas[c][d] =  kas2[c][d] + kas1[c][d] ;	
		
 		printf("%f\t", kas[c][d]);
 		
 }
 	 
	  printf("\n");
}


//******************************   RHS     ********************************


   printf("Enter the Heat flux in and out of the rod \n");
 
   for (c = 0; c <=e; c++){
   	
   	  for (d = 0; d < 1; d++){
   	  	
         scanf("%f", &q[c][d]);
		 }
   }

       for (c=0 ; c<=e; c++){
     	for (d =0; d<1;d++){
 		
 		printf("Thermal forces at node %d = %.3f\t",c+1, area*q[c][d]);
 }
  printf("\n");
  }
  

//****************************************Jacobi iterations**********************************
//Jacobi iteration to find coefficients form = kas*T = q ; We need to find T

printf("Enter the first approximation for temperature \n");
for(i=0;i<e+1;i++)
scanf("%f",&approxT[i][0]);

for(i=0;i <=e;i++)//We calculate the diagonal inverse matrix make all other entries as zero except Diagonal entries whose resciprocal we store
    for(j=0;j<=e ;j++)
    {    if(i==j)
       
       Diaginv[i][j]=(1./kas[i][j]);
        else
        Diaginv[i][j]=0;
        //printf("Diaginv [%d][%d] = %.2f\n",i,j, Diaginv[i][j]);
    }
    
    
for(i=0;i <=e;i++)
    for(j=0;j<=e ;j++)//calculating the E matrix L+U
    {     if(i==j)
        E[i][j]=0;
        else
        if(i!=j)
        E[i][j]=kas[i][j];
       // printf("L+U[%d][%d] = %f\n",i,j, E[i][j]);
    }

int iter;
printf("Enter the number of iterations:\n");
scanf("%d",&iter);//enter the number of iterations
int ctr=1;
int octr;
while(ctr<=iter)
{

	for(i=0; i<=e; i++)
	{
		for(j=0; j<1; j++)
		{
			
			sum=0;
			for(o=0; o<=e; o++)
			{
				sum = sum + E[i][o] * approxT[o][j];
			}
			F[i][j] = sum;
		}
	}
	
	for(i=0; i<=e; i++)
	{
		for(j=0; j<1; j++)
		{
		//	printf("(L+U)*ApproxT = [%d][%d] = %f\t",i,j,F[i][j]);
		}
		//printf("\n");
	}

      	for(i=0;i<=e;i++) //the matrix(RHS-Ex)
	// for(j=0;j<1;j++)
	 {
	 		RHS1[i][0] = q[i][0] - F[i][0] ;
	
//	printf("RHS-L+U*x = RHS1 = [%d][%d] = %f\n",i,j, RHS1[i][j]);
	 }
      
	for(i=0; i<=e; i++)
	{
		for(j=0; j<1; j++)
		{
				
			sum=0;
			for(o=0; o<=4; o++)
			{
					
				sum = sum + Diaginv[i][o] * RHS1[o][j];
			}
			RHS2[i][j] = sum;
		}
	}

	for(i=0; i<=e; i++)
	{
		for(j=0; j<1; j++)
		{
		//	printf("RHS2 = XN+1 = [%d][%d] = %.3f\t",i,j,RHS2[i][j]);
		}
		//printf("\n");
	}
      
      
for(octr=0;octr<=e;octr++)
    approxT[octr][0]= RHS2[octr][0];   //store x value in the next approximation
    printf("The Value after iteration %d is\n",ctr);
for(i=0;i<=e;i++)
for(j=0;j<1;j++)
printf("The Temperature at node %d = %.3f\n",i+1,approxT[i][j]);//display the value after the pass
ctr++;
}
		
	
getch();



return 0;
}	
	
