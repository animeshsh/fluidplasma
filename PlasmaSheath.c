//Animesh Sharma
// 
//MacCormack Method
//
// gcc -Wall -std=c99 -o projMac projMac.c -lm
//
// ./projMac>ProjMac.dat
//

#include<stdio.h>
#include<stdlib.h>
#include<math.h>

// Global variables-----------------------------------------------------------
const int rows=2;
const int npts=501; //101,1001
const int jmax=501; //101,1001
const int tmax=10001;
const int beta=1.8;
const double e4=-0.2;//0.2 0.4
int jj=0;

double kB=1.38064852e-23,e=1.60217662e-19,epo= 8.854187817e-12,pi=3.141549625;
double mi=2.18017e-25,Ti=300, Te=11600;
double po,To=300,no=1e14;

double *ui,*ni,*phi,*Ex,*temp;

double *U[2],*Unew[2],*F[2],*S[2];

double uB,lamD,tao,R,z,a;
double dx,dt;

/*
int r = 2, c = npts, i, j, count;
 
    int *arr[r];
    for (i=0; i<r; i++)
         arr[i] = (int *)malloc(c * sizeof(int));
	 
	for(int i = 0; i < N; i++)
		free(ptr[i]);
	free(ptr);
*/



// Function definitions---------------------------------------------------------------
void setval(double mat[],double val){
	for(int j=0;j<jmax;j++)
	{
		mat[j]=val;
	}
}

// set B[j]=A[j]*val
void setmatval(double B[],double A[],double val) {
	for(int j=0;j<jmax;j++)
	{
		B[j]=A[j]*val;
	}
}

// set C[j]=A[j]*B[j]*val
void setmatvalmul(double C[],double A[], double B[], double val) {
	for(int j=0;j<jmax;j++)
	{
		C[j]=A[j]*B[j]*val;
	}
}

// set C[j]=A[j]*B[j]*val
void setmatvaldiv(double C[],double A[], double B[], double val) {
	for(int j=0;j<jmax;j++)
	{
		C[j]=(A[j]/B[j])*val;
	}
}

void setmatvalsum(double C[],double A[], double B[]){	
	for(int j=0;j<jmax;j++)
	{
		C[j]=A[j]+B[j];
	}
}

void setup(){
	uB=sqrt(kB*Te/mi);
	lamD=pow(epo*kB*Te/(e*e*no),0.5);	
	R=100*lamD;
	dx=R/(double)(npts-1)/lamD;
	dt=0.01;//1e-8/(lamD/uB);
	a=0.;//a=0 , a=0.1;
	z=0.9/100;//0.85/100, 0.53/100
	tao=Ti/Te;
	printf("in setup \n");
    
	ui= malloc(jmax*sizeof(double));
	ni=malloc(jmax*sizeof(double));
	phi=malloc(jmax*sizeof(double));
    	Ex=malloc(jmax*sizeof(double));
	temp=malloc(jmax*sizeof(double));
	int r=rows;
	int c=jmax;
	for (int i=0; i<r; i++)
	{
         U[i] = (double *)malloc(c * sizeof(double));
		 Unew[i] = (double *)malloc(c * sizeof(double));
		 F[i] = (double *)malloc(c * sizeof(double));
		 S[i] = (double *)malloc(c * sizeof(double));	 
	}	  
}

void cleanup()
{
   
	free(ui);
	free(ni);
	free(phi);
    	free(Ex); 
	free(temp);
	int N=rows;
	for(int i = 0; i < N; i++)
	{
		free(U[i]);
		free(F[i]);
		free(S[i]);
		free(Unew[i]);
	}
}

void initialize(){
    	printf("in intialsize \n");
	setval(ui,0.);
	setval(phi,0.);
	setval(ni,1.);	
}

void tridiag(double *aa, double *dd, double *cc, double *bb, double *x, int imax) {
  //
  // Thomas's Tridiagonal Algorithm
  //
  // Description:
  //
  //   Solve a tridiagonal system of the form:
  //
  //   d[0] x[0] + c[0] x[1] = b[0]
  //   a[1] x[0] + d[1] x[1] + c[1] x[2] = b[1]
  //   a[2] x[1] + d[2] x[2] + c[2] x[3] = b[2]
  //   ...
  //   a[N-2] x[N-3] + d[N-2] x[N-2] + c[N-2] x[N-1] = b[N-2]
  //   a[N-1] x[N-2] + d[N-1] x[N-1] = b[N-1]
  //
  // Reference:
  //    Cheney and Kincaid (1994), Numerical Mathematics and Computing,
  //    3rd ed., Brooks/Cole Publishing Co., Pacific Grove CA, 1994,
  //    sec. 6.3, pp. 249-253
  //
  int i;
  double xmult;

  for (i = 1; i < imax; i++) {
    xmult = aa[i]/dd[i-1];
    dd[i] = dd[i] - xmult*cc[i-1];
    bb[i] = bb[i] - xmult*bb[i-1];
  }

  x[imax-1] = bb[imax-1]/dd[imax-1];

  for (i = imax-2; i >= 0; i--) {
    x[i] = (bb[i] - cc[i]*x[i+1])/dd[i];
  }
}


void calcphi(double *V, double *Ne){
	
	double *Vold = malloc(jmax*sizeof(double));
	double *aa = malloc(jmax*sizeof(double));
	double *bb = malloc(jmax*sizeof(double));
	double *cc = malloc(jmax*sizeof(double));
	double *dd = malloc(jmax*sizeof(double));
	setmatval(Vold,V,1);
	
	
	for(int k=0;k<10;k++)
	{
		Vold[0]=-50;
		Vold[jmax-1]=0;
		// unew[0] = uold[0]
		dd[0] = -(2.0+dx*dx*exp(Vold[0]));
		cc[0] = 1.;
		bb[0] =-Vold[0];

		// aa[i]*unew[i-1] + dd[i]*unew[i] + cc[i]*unew[i+1] = uold[i]


		for (int j = 1; j < jmax-1; j++) {
		aa[j] = 1.;
		dd[j] = -(2.0+dx*dx*exp(Vold[j]));
		cc[j] = 1.;
		bb[j] = -(dx*dx)*(ni[j]-exp(Vold[j]) + Vold[j]*exp(Vold[j]));
		}

		// unew[jmax-1] = uold[jmax-1]
		aa[jmax-1] = 1.;
		dd[jmax-1] = -(2.0+dx*dx*exp(Vold[jmax-1]));
		bb[jmax-1] =0.0;

		tridiag(aa, dd, cc, bb, V,jmax);
		for(int j=0;j<jmax;j++)
		{
			double t= V[j]-Vold[j];
			Vold[j]=(1-beta)*V[j]+beta*t;
		}
	}		
}

void calcEx(double *E, double *V){
	for(int j=0;j<=jmax-1;j++)
		E[j]=-(V[j+1]-V[j])/dx;	
	V[jmax-2]=0;
}

// main()--------------------------------------------------------------------------
int main(void)
{	
	FILE *fp;
	fp=fopen("proj2.dat","w");
	
	int j=0,n=0,i=0;
	
	printf("before setup\n");
	
	setup();
	
    printf("check1 \n");
    
	initialize();
	
    printf("check2 \n");

	calcphi(phi,ni);
	calcEx(Ex,phi);
    
    printf("check3 \n");
    
    for(j=0;j<jmax;j++)
    {
        U[0][j]=ni[j];
        U[1][j]=ni[j]*ui[j];
        F[0][j]=ni[j]*ui[j];
        F[1][j]=ni[j]*ui[j]*ui[j]+ni[j]*tao;
        S[0][j]=z*exp(phi[j]);
        S[1][j]=ni[j]*Ex[j]-ni[j]*fabs(ui[j])*ui[j]*a;
    }

    printf("check4 \n");
	
    fprintf(fp,"n:%d \n npts:%d\n alpha:%f\n t:%.4E\n uB:%.4E\n lamD:%.4E\n Ti/Te:%.4E\n,z:%.4E\n,dt:%.4E\n",0,npts,a,dt*(double)(tmax-1),uB,lamD,tao,z,dt);
	/*
	printf("%10s , %10s , %10s , %10s , %10s \n","j","x/lamDo","ni/no","ui/uB","phi");
	for(j=0;j<jmax;j++)
	{
		printf("%10d , %10.4E , %10.4E , %10.4E , %10.4E \n",j,(j*dx),ni[j],ui[j],phi[j]);
	}*/
    
	double *Ustar[2];
	double *Fstar[2];
	double *Sstar[2];
	double *phistar;
	double *Estar;
	double D4U;
	double nistar;
	double uistar;
	
	for (int i=0; i<rows; i++)
	{
         Ustar[i] = (double *)malloc(jmax * sizeof(double));
		 Fstar[i] = (double *)malloc(jmax * sizeof(double));
		 Sstar[i] = (double *)malloc(jmax * sizeof(double));	 
	}
	phistar = (double *)malloc(jmax * sizeof(double));
	Estar = (double *)malloc(jmax * sizeof(double));
	printf("check5\n");
	
	
	for(n=0 ; n < tmax ; n++)
	{
		if(n>2000);
			jj=3;
		calcphi(phi,ni);
		calcEx(Ex,phi);
		//ui[jmax-1]=z*ui[jmax-2]+dx;
		
		for(j=jj;j<jmax-jj;j++)
		{
			U[0][j]=ni[j];
			U[1][j]=ni[j]*ui[j];
			F[0][j]=ni[j]*ui[j];
			F[1][j]=ni[j]*ui[j]*ui[j]+ni[j]*tao;
			S[0][j]=z*exp(phi[j]);
			S[1][j]=ni[j]*Ex[j]-ni[j]*fabs(ui[j])*ui[j]*a;
		}
		for(j=jj;j<jmax-jj;j++)
		{
			double delU;
			
			//printf("check6\n");
			for(i=0;i<rows;i++)
			{
				if(j<=1) 			D4U=3.*U[i][j]-14.*U[i][j+1]+26.*U[i][j+2]-24.*U[i][j+3]+11.*U[i][j+4]-2.*U[i][j+5];
				if(j>1&&j<jmax-2) 	D4U=U[i][j-2]-4.*U[i][j-1]+6.*U[i][j]-4.*U[i][j+1]+U[i][j+2];
				if(j>=jmax-2) 		D4U=-2*U[i][j-5]+11.*U[i][j-4]-24.*U[i][j-3]+26.*U[i][j-2]-14.*U[i][j-1]+3.*U[i][j];
				
				delU=-(dt/dx)*(F[i][j+1]-F[i][j])+S[i][j]*dt+e4*D4U*dt;
				Ustar[i][j]=U[i][j]+delU;
			}
			nistar=Ustar[0][j];
			uistar=Ustar[1][j]/Ustar[0][j];
			Fstar[0][j]=nistar*uistar;
			Fstar[1][j]=nistar*uistar*uistar+nistar*tao;
	
		}
		//printf("check7\n");
		setmatval(phistar,phi,1);
		calcphi(phistar,Ustar[0]);
		calcEx(Estar,phistar);
		for(j=jj;j<jmax-jj;j++){
			nistar=U[0][j];
			uistar=U[1][j]/U[0][j];
			Sstar[0][j]=z*exp(phistar[j]);
			Sstar[1][j]=nistar*Estar[j]-nistar*fabs(uistar)*uistar*a;
		}
		for(j=jj;j<jmax-jj;j++)
		{
			for(i=0;i<rows;i++)
			{
				if(j<=1) 			D4U=3.*Ustar[i][j]-14.*Ustar[i][j+1]+26.*Ustar[i][j+2]-24.*Ustar[i][j+3]+11.*Ustar[i][j+4]-2.*Ustar[i][j+5];
				if(j>1&&j<jmax-2) 	D4U=Ustar[i][j-2]-4.*Ustar[i][j-1]+6.*Ustar[i][j]-4.*Ustar[i][j+1]+Ustar[i][j+2];
				if(j>=jmax-2) 		D4U=-2*Ustar[i][j-5]+11.*Ustar[i][j-4]-24.*Ustar[i][j-3]+26.*Ustar[i][j-2]-14.*Ustar[i][j-1]+3.*Ustar[i][j];
				
				double delUstar=-(dt/dx)*(Fstar[i][j]-Fstar[i][j-1])+Sstar[i][j]*dt+e4*D4U*dt;
				Unew[i][j]=0.5*(U[i][j]+Ustar[i][j]+delUstar);
			}
		}
		for(i=0;i<rows;i++)
			setmatval(U[i],Unew[i],1);
	
		setmatval(ni,U[0],1);
		setmatvaldiv(ui,U[1],U[0],1);
	
		/*
		printf("n:%d , t:%.4E , uB: %.4E,lamD: %.4E, Ti/Te:%.4E , z:%.4E \n",n,dt*(double)n,uB,lamD,tao,z);
		printf("%10s , %10s , %10s , %10s , %10s \n","j","x/lamDe","ni/no","ui/uB","phi");
		for(j=0;j<jmax;j++)
		{
			printf("%10d , %10.4E , %10.4E , %10.4E , %10.4E \n",j,(j*dx),ni[j],ui[j],phi[j]);
		}*/
	//	calcphi(phi,ni);
	//	calcEx(Ex,phi);
	//	for(j=1 ; j<5 ; j++);
	//		ni[j]=(ni[j-1]+ni[j+1])/2;
			
	//	for(j=jmax-6 ; j<jmax-1 ; j++)
	//		ni[j]=(ni[j-1]+ni[j+1])/2;
    }
	fprintf(fp,"%10s , %10s , %10s , %10s , %10s , %10s \n","j","x/lamDe","ni/no","ne/no","ui/uB","phi");
	for(j=0;j<jmax;j++)
		{
			fprintf(fp,"%10d , %10.4E , %10.4E , %10.4E , %10.4E, %10.4E \n",j,(j*dx),ni[j],exp(phi[j]),ui[j],phi[j]);
			printf("%10d , %10.4E , %10.4E , %10.4E , %10.4E, %10.4E \n",j,(j*dx),ni[j],exp(phi[j]),ui[j],phi[j]);
		}
	for(int i = 0; i < rows; i++)
	{
		free(Ustar[i]);
		free(Fstar[i]);
		free(Sstar[i]);
	}
   
	cleanup();
	fclose(fp);
	return 0;
}