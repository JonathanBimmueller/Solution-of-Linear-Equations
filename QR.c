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



void QRdecomposition(int n, int p, double* X, double *Q, double* R){
	for(int i=0; i<n; i++){
		for(int j=0; j<p; j++){
			if(i<p){				
				R[i*p+j]=0;
			}
			Q[i*p+j]=X[i*p+j];
		}
	}
	for(int j=0; j<p; j++){
		double x_j=Q[0*p+j]*Q[0*p+j];	//int
		for(int i=1; i<n; i++){
			x_j+=Q[i*p+j]*Q[i*p+j];
		}
		R[j*p+j]=sqrt(x_j);
		for(int i=0; i<n; i++){
			Q[i*p+j]=Q[i*p+j]/R[j*p+j];
		}
		for(int k=j+1; k<p; k++){
			R[j*p+k]=Q[0*p+j]*Q[0*p+k];
			for(int i=1; i<n; i++){
				R[j*p+k]+=Q[i*p+j]*Q[i*p+k];
			}
			for(int i=0; i<n; i++){		
				Q[i*p+k]-=R[j*p+k]*Q[i*p+j];	
			}			
		}
	}
}



void QRsolve(int n, int p, double* Q, double* R, double* y, double* x){
	for(int i=0; i<p; i++){
		x[i]=Q[0*p+i]*y[0];
		for(int j=1; j<n; j++){
			x[i]+=Q[j*p+i]*y[j];
		}
	}
	for(int i=(p-1);i>=0;i--){
		for(int j=(i+1);j<p;j++){
			x[i]-=R[i*p+j]*x[j];
		}
		x[i]=x[i]/R[i*p+i];
	}
}



int main(){
	int n=5;
	int p=3;
	
	double* X = malloc(n*p*sizeof(*X));
	X[0*p+0]=1;		X[0*p+1]=-2;		X[0*p+2]=1;
	X[1*p+0]=1;		X[1*p+1]=3;			X[1*p+2]=-1;
	X[2*p+0]=1;		X[2*p+1]=0;			X[2*p+2]=-1;
	X[3*p+0]=1;		X[3*p+1]=3;			X[3*p+2]=-1;
	X[4*p+0]=0;		X[4*p+1]=sqrt(7);	X[4*p+2]=11/sqrt(7);
	printf("X:\n");
	PrintMatrix(n,p,X);
	
	double* Q = malloc(n*p*sizeof(*X));
	double* R = malloc(p*p*sizeof(*X));
	
	QRdecomposition(n, p, X, Q, R);
	printf("QR-Zerlegung von X:\nQ:\n");
	PrintMatrix(n,p,Q);
	printf("R:\n");
	PrintMatrix(p,p,R);
	
	free(X);
	free(Q);
	free(R);
	
	n=3;
	p=3;
	
	X = malloc(n*p*sizeof(*X));
	X[0*p+0]=2;		X[0*p+1]=1;		X[0*p+2]=1;
	X[1*p+0]=4;		X[1*p+1]=2;		X[1*p+2]=-1;
	X[2*p+0]=-1;	X[2*p+1]=0;		X[2*p+2]=7;
	printf("X:\n");
	PrintMatrix(n,p,X);
	
	Q = malloc(n*p*sizeof(*X));
	R = malloc(p*p*sizeof(*X));
	
	QRdecomposition(n, p, X, Q, R);
	printf("QR-Zerlegung von X:\nQ:\n");
	PrintMatrix(n,p,Q);
	printf("R:\n");
	PrintMatrix(p,p,R);
	
	double* y = malloc(n*sizeof(*y));
	y[0]=1;	y[1]=2;	y[2]=3;
	
	double* b = malloc(p*sizeof(*b));
	
	QRsolve(n,p,Q,R,y,b);
	printf("Loesung des gleichen Gleichungsystems wie in der LU-Zerlegung:\n");
	PrintVector(p,b);
	
	free(X);
	free(Q);
	free(R);
	free(y);
	free(b);
	
	n=6;
	p=2;
	
	X = malloc(n*p*sizeof(*X));
	X[0*p+0]=0.5;	X[0*p+1]=1;
	X[1*p+0]=1.5;	X[1*p+1]=1;
	X[2*p+0]=3.5;	X[2*p+1]=1;
	X[3*p+0]=5;		X[3*p+1]=1;
	X[4*p+0]=7.5;	X[4*p+1]=1;
	X[5*p+0]=8;		X[5*p+1]=1;
	printf("X:\n");
	PrintMatrix(n,p,X);	

	Q = malloc(n*p*sizeof(*X));
	R = malloc(p*p*sizeof(*X));
	
	QRdecomposition(n, p, X, Q, R);
	printf("QR-Zerlegung von X:\nQ:\n");
	PrintMatrix(n,p,Q);
	printf("R:\n");
	PrintMatrix(p,p,R);
	
	y = malloc(n*sizeof(*y));
	y[0]=160;	y[1]=140;	y[2]=120;	y[3]=120;	y[4]=100;	y[5]=80;	y[5]=80;
	
	b = malloc(p*sizeof(*b));
	
	QRsolve(n,p,Q,R,y,b);
	printf("Loesung des linearen Ausgleichproblems aus Numerik I:\n");
	PrintVector(p,b);
	printf("Dies bedeutet:\n f(x)=%lf*x+%lf\n mit Nullstelle %lf\n\nDer Professor glaubt also seine Vorlesung nach %d Wochen leergelesen zu haben",b[0],b[1],-b[1]/b[0],(int)(-b[1]/b[0]));
	
	
	free(X);
	free(Q);
	free(R);
	free(y);
	free(b);
}