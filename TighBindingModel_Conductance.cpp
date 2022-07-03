# include <iostream>
# include <math.h>
#include "mkl_lapack.h"
#include "mkl_blas.h"
const int SIZE= 12;
const int M=5;
using namespace std;

main()
{
int i, j, k, l, dif, r, N;
N=SIZE*M;
double E=-0.5, EF=0.0, V=0.0, t=1.0, t1=0.1, t3=0.1;

char trans='T';
char TRANS='N';
char transa='N';
char transb='N';
char transH='T';

double alpha=1.0;
double beta=0.0;
double trace=0.0;

double G[N][N];
double GT[N][N];
double C[N][N];
double Q[N][N];
double P[N][N];
double gamaL[N][N];
double gamaR[N][N];
double work[N][N];

double HC[SIZE][SIZE];
double a[SIZE][SIZE];
double A[SIZE][SIZE];
double b[SIZE][SIZE];
double H00[SIZE][SIZE];
double H10[SIZE][SIZE];
double H01[SIZE][SIZE];
double B[SIZE][SIZE];
double D[SIZE][SIZE]={0.0};
double d[SIZE][SIZE]={0.0};
double c[SIZE][SIZE]={0.0};
double T[SIZE][SIZE]={0.0};
double TB[SIZE][SIZE]={0.0};
double L[SIZE][SIZE]={0.0};
double R[SIZE][SIZE]={0.0};



MKL_INT n=SIZE;
MKL_INT ipiv[SIZE];
MKL_INT info;
MKL_INT lda=SIZE;
MKL_INT ldb=SIZE;
MKL_INT ldc=SIZE;
MKL_INT nrhs = SIZE;

MKL_INT lwork=N;
MKL_INT IPIV[N];
MKL_INT LDA=N;



/***************************  E-H(00)  *****************************/
for(i=0;i<SIZE;i++){
	for(j=0;j<=i;j++){
		dif=i-j;
		switch(dif){
		case(0):
		if((i+1)%4==1 || (i+1)%4==2)H00[i][i]=E;
		else H00[i][i]=E;break;
		case(1):
			 r=(j+1)%4;
			 switch (r){
				 case 0:
					H00[i][j]=t3;
					H00[j][i]=t3;break;
				 case 2:
					H00[i][j]=t1;
					H00[j][i]=t1;break;
				 case 1:
					H00[i][j]=t;
					H00[j][i]=t;break;
				 default:
					H00[i][j]=t;
					H00[j][i]=t;
					}break;
		case(3):
			r=(j+1)%4;
			switch (r){
				 case 0:
					H00[i][j]=t;
					H00[j][i]=t;
				 case 2:
					H00[i][j]=t;
					H00[j][i]=t;break;
				 case 1:
					H00[i][j]=t3;
					H00[j][i]=t3;break;
				 default:
					H00[i][j]=0.0;
					H00[j][i]=0.0;
					}break;
		default :
			 H00[i][j]=0.0;
			 H00[j][i]=0.0;
			}
			}
		}
cout<<"************  Matrix E-H(00) is: ************"<<"\n\n";
for(i=0;i<SIZE;i++){
	cout<<"\n\n";
	for(j=0;j<SIZE;j++){
		A[i][j]=H00[i][j];
      		a[i][j]=H00[i][j];
		cout<<A[i][j]<<"\t";
				}
			}
cout<<"\n\n\n";

/*************************** H(01) *****************************/

for (i=0;i<SIZE;i++){
	for(j=0;j<SIZE;j++){
		dif=i-j;
		switch(dif){
			case(0):
				H01[i][i]=0.0;break;

			case(1):
				r=(j+1)%8;
				switch (r){
					case 1:
						H01[i][j]=t;
						break;
					case 7:
						H01[i][j]=t;
						break;
					case 4:
						H01[i][j]=t3;
						break;
					default:
						H01[i][j]=0.0;

						}break;
			case(-1):

				r=(j+1)%8;
				switch (r){
					case 4:
						H01[i][j]=t;

					case 6:
						H01[i][j]=t;
						break;
					case 1:
						H01[i][j]=t3;
						break;
					default:
						H01[i][j]=0.0;
						}break;
			default :
				H01[i][j]=0.0;


				}
			}
		}
 cout<<"************   H(01) is:  ************"<<"\n\n";
for (i=0;i<SIZE;i++){
	for(j=0;j<=SIZE-1;j++){
			B[i][j]=H01[i][j];
 			cout<<B[i][j]<<"\t";
				}
 	cout<<"\n\n";
			}
cout<<"\n\n\n";

/********************** Complex Conjugate of H(01) ******************/

//cout<<"********* Complex Conjugate of H(01) is:  **********"<<"\n\n";
for (i=0;i<SIZE;i++){
	for(j=0;j<SIZE;j++){
			b[i][j]=B[j][i];
			H10[i][j]=b[i][j];
			//cout<<b[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}
cout<<"\n\n";

/************************  Calculation of t(0)=B  *********************/

dgetrf( &n, &n, &a[0][0], &lda, ipiv, &info ) ;
dgetrs( &trans, &n, &nrhs, &a[0][0], &lda, ipiv, &B[0][0], &ldb, &info );

//Check for success
if(info == 0){
	// Write the answer
	cout<<"************** t(0) is: **************"<<"\n\n";
	for (i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			T[i][j]=B[i][j];
			D[i][j]=B[i][j];
			cout<<B[j][i]<<"\t";
			
 				}
		cout<<"\n\n";
			}
	}
else {
	// Write an error message
	cout << "dgetrs returned error " << info << "\n";
	}

		// return info;


/***************************  Calculation of  t(0)tilda=b  ***********************/

for(i=0;i<SIZE;i++){
	for(j=0;j<SIZE;j++){
		a[i][j]=A[i][j];
				}
			}

dgetrf( &n, &n, &a[0][0], &lda, ipiv, &info ) ;
dgetrs( &trans, &n, &nrhs, &a[0][0], &lda, ipiv, &b[0][0], &ldb, &info );

//Check for success
		  if(info == 0){
			// Write the answer
			cout<<"***********t(0)tilda is:***********"<<"\n\n";
			for (i=0;i<SIZE;i++){
				for(j=0;j<SIZE;j++){
					d[i][j]=b[i][j];
					TB[i][j]=b[i][j];
					cout<<b[j][i]<<"\t";
					
						 }
				cout<<"\n\n";
					}
					 }
		  else{
			// Write an error message
			cout << "dgetrs returned error " << info << "\n";
			 }

			// return info;



for(l=1;l<=12;l++){
//cout<<"\n\n\n";
//A=a=1-t(l-1)*ttilda(l-1)-ttilda(l-1)*t(l-1)
//cout<<l<<"\t\t1-t("<<l-1<<")t("<<l-1<<")tilda-t("<<l-1<<")tildat("<<l-1<<") is  :"<<"\n\n\n";
	for (i=0;i<SIZE;i++){
		for (j=0;j<SIZE;j++){
			A[i][j]=0.0;
			a[i][j]=0.0;
		for (k=0;k<SIZE;k++){
			A[i][j]+=-B[i][k]*b[k][j]-b[i][k]*B[k][j];
				}
			if(i==j)A[i][i]=A[i][i]+1.0;
			a[i][j]=A[i][j];
		//cout<<a[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}

dgemm(&transa, &transb, &n, &n, &n, &alpha, &B[0][0], &lda, &B[0][0], &ldb, &beta, &c[0][0], &ldc);//t(l-1)^2=c[][]

dgetrf( &n, &n, &a[0][0], &lda, ipiv, &info ) ;
dgetrs( &TRANS, &n, &nrhs, &a[0][0], &lda, ipiv, &c[0][0], &ldb, &info ); //t(l)=c[][]
//Check for success
if(info == 0){
	// Write the answer
	//cout<<"*********  t("<<l<<") is: *******"<<"\n\n";
	for (i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			B[i][j]=c[i][j];		//t(l)=B[][]
			//cout<<B[i][j]<<"\t";
				 }
		//cout<<"\n\n";
				}
		}
 else{
	// Write an error message
	cout << "dgetrs returned error " << info << "\n";
	}
	// return info;

dgemm(&transa, &transb, &n, &n, &n, &alpha, &d[0][0], &lda, &B[0][0], &ldb, &beta, &c[0][0], &ldc);//t(l)*ttilda(l-1)=c[][]
//cout<<"****************T("<<l<<")***************\n";
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		T[i][j]+=c[i][j];
		//cout<<T[i][j]<<"\t";
				}
  	//cout<<"\n\n";
			}
dgemm(&transa, &transb, &n, &n, &n, &alpha, &b[0][0], &lda, &b[0][0], &ldb, &beta, &c[0][0], &ldc);//ttilda(l-1)^2=c[][]
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
			a[i][j]=A[i][j];
				}
		}
dgetrf( &n, &n, &a[0][0], &lda, ipiv, &info ) ;
dgetrs( &TRANS, &n, &nrhs, &a[0][0], &lda, ipiv, &c[0][0], &ldb, &info ); //ttilda(l)=b[][]

	//Check for success
	if(info == 0){
	// Write the answer
	//cout<<"       t("<<l<<")tilda is:      "<<"\n\n";
	for (i=0;i<SIZE;i++){
		for(j=0;j<SIZE;j++){
			b[i][j]=c[i][j];
			//cout<<b[i][j]<<"\t";
				 }
		//cout<<"\n\n";
				}
		 }
	else{
	// Write an error message
	cout << "dgetrs returned error " << info << "\n";
	}

	// return info;

dgemm(&transa, &transb, &n, &n, &n, &alpha, &D[0][0], &lda, &b[0][0], &ldb, &beta, &c[0][0], &ldc);//tilda(l)*t(l-1)=c[][]
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		TB[i][j]+=c[i][j];
	//	cout<<TB[i][j]<<"\t"; //TB
				}
  //	cout<<"\n\n";
			}

dgemm(&transa, &transb, &n, &n, &n, &alpha, &D[0][0], &lda, &B[0][0], &ldb, &beta, &c[0][0], &ldc);//t(l)*t(l-1)=c[][]
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		D[i][j]=0.0;
		D[i][j]=c[i][j];
	//	cout<<d[i][j]<<"\t"; //t(l)*t(l-1)=C[][]=d[][]
				}
  //	cout<<"\n\n";
			}

dgemm(&transa, &transb, &n, &n, &n, &alpha, &d[0][0], &lda, &b[0][0], &ldb, &beta, &c[0][0], &ldc);//ttilda(l-1)*ttilda(l)=c[][]
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		d[i][j]=0.0;
		d[i][j]=c[i][j];
	//	cout<<d[i][j]<<"\t"; //ttilda(l)*tt(l-1)=C[][]=d[][]
				}
  //	cout<<"\n\n";
			}
			}//end of l

// cout<<"       ******************* T is: *******************"<<"\n\n";
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
// 		cout<<T[i][j]<<"\t";
				}
// 	cout<<"\n\n";
			}
// cout<<"\n\n";

// cout<<"       ******************* TB is: *******************"<<"\n\n";
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
// 		cout<<TB[i][j]<<"\t";
				}
// 	cout<<"\n\n";
			}
// cout<<"\n\n";

dgemm(&transH, &transb, &n, &n, &n, &alpha, &H01[0][0], &lda, &T[0][0], &ldb, &beta, &c[0][0], &ldc);//H01*T=c[][]
// cout<<"       ******************* sigma R is: *******************"<<"\n\n";
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		R[i][j]=c[i][j];
// 		cout<<R[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}
//cout<<"\n\n";
//cout<<"       ******************* gama R is: *******************"<<"\n\n";
for (i=0;i<N;i++){
	for (j=0;j<N;j++){
		if(i>((M-1)*SIZE) && j>((M-1)*SIZE))gamaR[i][j]=R[i-(M-1)*SIZE][j-(M-1)*SIZE]-R[j-(M-1)*SIZE][i-(M-1)*SIZE];
		else gamaR[i][j]=0.0;
		//cout<<gamaR[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}
//cout<<"\n\n";


dgemm(&transH, &transb, &n, &n, &n, &alpha, &H10[0][0], &lda, &TB[0][0], &ldb, &beta, &c[0][0],&ldc);//H10*TB=c[][]
//cout<<"       ******************* sigma L is: *******************"<<"\n\n";
for (i=0;i<SIZE;i++){
	for (j=0;j<SIZE;j++){
		L[i][j]=c[i][j];
		//cout<<L[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}
//cout<<"\n\n";
//cout<<"       ******************* gama L is: *******************"<<"\n\n";
for (i=0;i<N;i++){
	for (j=0;j<N;j++){
		if(i<SIZE && j<SIZE)gamaL[i][j]=L[i][j]-L[j][i];
		else gamaL[i][j]=0.0;
		//cout<<gamaL[i][j]<<"\t";
				}
	//cout<<"\n\n";
			}
//cout<<"\n\n";

/**********************************************************************************************************/
/**********************************************************************************************************/

for(EF=-0.04;EF<=0.04;EF+=0.001){
trace=0.0;
for (i=0;i<N;i++){
	for(j=0;j<N;j++){
		G[i][j]=0.0;
			}
		}
for(i=0;i<SIZE;i++){
	//cout<<"\n\n";
	for(j=0;j<SIZE;j++){
		if((i+1)%4==1 || (i+1)%4==2)H00[i][i]=EF+V;
		else H00[i][i]=EF-V;
		//cout<<H00[i][j]<<"\t";
				}
			}
//cout<<"\n\n";
			
		
for (k=0;k<M;k++){
	for (i=SIZE*k;i<SIZE*(k+1);i++){
		for(j=SIZE*k;j<SIZE*(k+1);j++){
			G[i][j]=H00[i-SIZE*k][j-SIZE*k];
			if(k==0)G[i][j]=H00[i-SIZE*k][j-SIZE*k];//-L[i-SIZE*k][j-SIZE*k];
			if(k==M-1)G[i][j]=H00[i-SIZE*k][j-SIZE*k];//-R[i-SIZE*k][j-SIZE*k];
							}
					}
		}
for (k=0;k<M-1;k++){
	for (i=SIZE*k;i<SIZE*(k+1);i++){
		for(j=SIZE*(k+1);j<SIZE*(k+2);j++){
			G[i][j]=-H01[(i-SIZE*k)][(j-SIZE*(k+1))];
			G[j][i]=-H10[(i-SIZE*k)][(j-SIZE*(k+1))];
						}
					}
		}

//cout<<"E="<<E<<"\t\t    ************   E - H center- sigma L -sigma R is:  ************"<<"\n\n";
for (i=0;i<N;i++){
	for(j=0;j<N;j++){
//			cout<<G[i][j]<<"\t";
				}
//	cout<<"\n\n";
			}
//cout<<"\n";
		

dgetrf( &N, &N, &G[0][0], &LDA, IPIV, &info ) ;
dgetri( &N, &G[0][0], &LDA, IPIV, &work[0][0], &lwork, &info );//G=inverse G
//Check for success
	if(info == 0){
		// Write the answer
		//cout<<"      ***************  G is:   ***************   "<<"\n\n";
		for (i=0;i<N;i++){
			for(j=0;j<N;j++){
				GT[i][j]=G[j][i];
				//cout<<G[i][j]<<"\t";
					}
			//cout<<"\n\n";
				}
		}
	else {
		// Write an error message
		cout << "G has not invers " << info << "\n";
		}

		// return info;

dgemm(&transa, &transb, &N, &N, &N, &alpha, &gamaL[0][0], &LDA, &G[0][0], &LDA, &beta, &C[0][0],&LDA);//gamaL*G=C[][]
dgemm(&transa, &transb, &N, &N, &N, &alpha, &gamaR[0][0], &LDA, &GT[0][0], &LDA, &beta, &P[0][0],&LDA);//GamaR*GT=G[][]
dgemm(&transa, &transb, &N, &N, &N, &alpha, &C[0][0], &LDA, &P[0][0], &LDA, &beta, &Q[0][0],&LDA);//gamaL*G*gamaR*G(transpose)=GT[][]

for(i=0;i<N;i++){
	//for(j=0;j<N;j++){
		//cout<<Q[i][j]<<"\t";
		trace+=Q[i][i];
		}
	//cout<<"\n\n";
	//	}		
//cout<<EF<<"\t\t"<<trace<<"\n";
		}


return 0;
}
