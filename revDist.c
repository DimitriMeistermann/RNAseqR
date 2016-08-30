#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

//Algorithmes de calcul de distance de chaîne de caractère

SEXP  orderDist(SEXP vA, SEXP vB){
    //input : two integer vector
	int n;
    int i = 0;
	int* a;
	int* b;
	long res = 0;
	int* compTab;
    SEXP resR;
	
	n=length(vA);
    
	a = INTEGER(vA);
	b = INTEGER(vB);
	
	compTab=Calloc(n,int);
	PROTECT(resR=allocVector(INTSXP,1));
    
	for(i=0;i<n;i++) compTab[b[i]-1]=i;
	for(i=0;i<n;i++) res+= abs(compTab[a[i]-1]-i);
	
	INTEGER(resR)[0]=res;
	Free(compTab);
	UNPROTECT(1);
	return resR;
}


SEXP  LevDist(SEXP vA, SEXP vB){
    //LevensteinDistance
	//input : two integer vector
	int ncols;
	int nrows;
    int i = 0;
	int j = 0;
	int* a;
	int* b;
	long** Mat;
	long del;
	long ins;
	long sub;
	long cost;
    SEXP resR;
	
	nrows=length(vA);
	ncols=length(vB);
	
	a = INTEGER(vA);
	b = INTEGER(vB);

	Mat=Calloc(nrows,long*);
	for(i=0;i<=nrows;i++){
		*Mat = Calloc(ncols,long);
		Mat[i][0]=i;
	}
	for(j=0;j<=ncols;j++) Mat[0][j]=j;

	PROTECT(resR=allocVector(INTSXP,1));
    printf("ready2");
	for(i=1;i<=nrows;i++){
		for(j=1;j<=ncols;j++){
			cost = (a[i-1]==b[j-1]) ? 0 : 1;
			del=Mat[i-1][j]+1;
			ins=Mat[i][j-1]+1;
			sub=Mat[i-1][j-1]+cost;
			if(del>=ins && del>=sub) Mat[i][j]=del;
			if(ins>=del && ins>=sub) Mat[i][j]=ins;
			if(sub>=ins && sub>=del) Mat[i][j]=sub;
		}
	}
	INTEGER(resR)[0]=Mat[nrows][ncols];
	
	for(i=0;i<=nrows;i++) Free(Mat[i]);
	Free(Mat);
	UNPROTECT(1);
	return resR;
}
