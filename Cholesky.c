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



void Choleskydecomposition (int n, double* A){
	for(int k=0; k<n; k++){
		A[k*n+k]=sqrt(A[k*n+k]);
		for(int j=k+1; j<n; j++){
			A[j*n+k]=A[j*n+k]/A[k*n+k];
		}
		for(int i=k+1; i<n; i++){
			for(int j=i; j<n; j++){
				A[j*n+i]-=A[i*n+k]*A[j*n+k];
			}
		}
	}
}



int main(){
	int n=3;
	int p=3;
	
	double* A = malloc(n*p*sizeof(*A));
	A[0*n+0]=1;		A[0*n+1]=2;		A[0*n+2]=3;
	A[1*n+0]=2;		A[1*n+1]=13;	A[1*n+2]=18;
	A[2*n+0]=3;		A[2*n+1]=18;	A[2*n+2]=26;
	printf("A:\n");
	PrintMatrix(n,n,A);
	
	Choleskydecomposition(n, A);
	printf("Cholesky-Zerlegung von A:\n");
	PrintMatrix(n,n,A);
	
	free(A);
}