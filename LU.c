#include <stdio.h>
#include <stdlib.h>
#include <math.h>



void PrintMatrix(int n, int p, double* X){
	for(int i=0; i<n; i++){
		for(int j=0; j<p; j++){
			printf("%lf ", X[i*p+j]);
		}
		printf("\n");
	}
	printf("\n");
}

void PrintVector(int n, double* x){
	for(int i=0; i<n; i++){
		printf("%lf \n", x[i]);
	}
	printf("\n");
}



void LUdecomposition(int n, double* A, double* perm){
	for(int k=0; k<n; k++){
		perm[k]=(double)k;
	}
	// find pivot row
	for(int k=0; k<n; k++){
		int p=k;
		for(int i=k+1; i<n; i++){
			if(fabs(A[i*n+k])>fabs(A[p*n+k])){
				p=i;
			}
		}
		// interchange permutation and rows
		double x=perm[k];
		perm[k]=perm[p];
		perm[p]=x;
		for(int j=0; j<n; j++){
			double y=A[p*n+j];
			A[p*n+j]=A[k*n+j];
			A[k*n+j]=y;
		}
		// compute multipliers
		for(int i=k+1; i<n; i++){
			A[i*n+k]=A[i*n+k]/A[k*n+k];
		}
		// adjust multipliers
		for(int i=k+1; i<n; i++){
			for(int j=k+1; j<n; j++){
				A[i*n+j] -= A[i*n+k]*A[k*n+j];
			}
		}
	}
}



void LUsolve(int n, double* A, double* perm, double* b, double* z, double* x){
	for(int i=0; i<n; i++){
		z[i]=b[(int)perm[i]];
		for(int j=0; j<=i-1; j++){
			z[i] -= A[i*n+j]*z[j];
		}
	}
	for(int i=n-1; i>=0; i--){
		x[i]=z[i];
		for(int j=i+1; j<n; j++){
			x[i] -= A[i*n+j]*x[j];
		}
		x[i]=x[i]/A[i*n+i];
	}
}



int main(){
	int n=3;
	int p=3;
	
	double* A = malloc(n*p*sizeof(*A));
	A[0*n+0]=2;		A[0*n+1]=1;		A[0*n+2]=1;
	A[1*n+0]=4;		A[1*n+1]=2;		A[1*n+2]=-1;
	A[2*n+0]=-1;	A[2*n+1]=0;		A[2*n+2]=7;
	printf("A:\n");
	PrintMatrix(n,p,A);
	
	double* perm = malloc(n*sizeof(*perm));
	
	LUdecomposition(n, A, perm);
	printf("LU-Zerlegung von A:\n");
	PrintMatrix(n,p,A);
	printf("perm:\n");
	PrintVector(n, perm);
	
	double* b = malloc(n*sizeof(*b));
	b[0]=1;	b[1]=2;	b[2]=3;
	printf("b:\n");
	PrintVector(p,b);
	
	double* z = calloc(n,sizeof(*z));
	double* x = calloc(n,sizeof(*x));
	
	LUsolve(n, A, perm, b, z, x);
	printf("Loesung x des Gleichungssystem Ax=b:\n");
	PrintVector(p,x);
	
	free(A);
	free(b);
	free(perm);
	free(z);
	free(x);
}