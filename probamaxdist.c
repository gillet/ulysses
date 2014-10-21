#include <time.h>
#include <R.h>
#include <Rmath.h>

//Simple matrix allocation
double **CallocMat(int nrow, int ncol){
  int i ;
  double *bigvec = Calloc( nrow*ncol , double) ;
  double **mat = Calloc( nrow, double *) ;
  for (i = 0 ; i < nrow ; i++){
    mat[i] = bigvec + i * ncol ;
  }
  return(mat) ;
}

//reset all values of the matrix to 0 
void ClearMat(double **mat, int nrow, int ncol){
  int i,j ;
  for (int i = 0 ; i < nrow ; ++i){
    for (int j = 0 ; j < ncol ; ++j){
      mat[i][j] = 0.0 ;
    }
  }
}


double MatSum(double **mat, int nrow, int ncol){
  int i,j;
  double sum = 0.0 ;
  for (i = 0 ; i < nrow; ++i)
    for(j= 0 ; j < ncol ; ++j)
      sum += mat[i][j];
  return(sum) ; 
}

void AfficheMat(double **mat, int nrow, int ncol){
  int i,j ;
  for (i = 0 ; i < nrow; ++i)
    {
      for(j= 0 ; j < ncol ; ++j)
	Rprintf("%.e ", mat[i][j]);
      Rprintf("\n") ;
    }
}

/* Proba of 2 consistent elements at max distance IScn given the distribution over distances
   ISdist */ 
void* Proba2Consistent(double *ISdist, double** pmat, int ndist, int IScn){
  int i,j,jmax ;
  for (i = 0; i < ndist ; ++i){
    pmat[i][i] = ISdist[i] * ISdist[i] ;
    jmax = imin2(ndist, i + IScn)  ;
    for (j = i+1; j <  jmax ; ++j){
      pmat[i][j] = 2 * ISdist[i] * ISdist[j];
    }
  }
}



/* Proba of the next 
*/
void ProbaNextConsistent(double *ISdist, double** pmat_next, double** pmat, int ndist, int IScn)
{
  int i, j, k, jmax ;
  double sij_j, si_ij, s_IS;
  for (i = 0 ; i < ndist ; ++i)
    {
    jmax = imin2(ndist, i + IScn);
    si_ij = pmat[i][i] ;
    s_IS = 0 ;
    
    for (j = i; j < jmax ; ++j )
      {
	sij_j = 0;
	if (j > i +1)
	  s_IS += ISdist[j-1] ;
	for (k = i ; k <= j; ++k) 
	  sij_j += pmat[k][j] ;
	si_ij += pmat[i][j];
	pmat_next[i][j] = ISdist[i] * sij_j + ISdist[j] * si_ij + pmat[i][j] * s_IS ;
      }
    }
}

void nPsdist(double* ISdist, double* probas, int *ndist, int *nPS, int *IScn){
  int i,j, n ;
  double** pmat; 
  double** pmat_next;
  //Rprintf("\n*******\nParameters, ndist: %d, nPS %d , IScn : %d,Allocating memory 1\n", *ndist, *nPS, *IScn) ;
  //probas = Calloc(*nPS, double) ;
  //Rprintf("Dynamic Programming for pIS up to n = %d\n", *nPS) ;
  probas[0] = 1; 
  pmat = CallocMat(*ndist, *ndist) ;    
  ClearMat(pmat, *ndist, *ndist) ;
  pmat_next = CallocMat(*ndist, *ndist) ;
  ClearMat(pmat_next, *ndist, *ndist) ;
  Proba2Consistent(ISdist, pmat, *ndist, *IScn) ;
  //Rprintf("#") ;
  probas[1] = MatSum(pmat, *ndist, *ndist) ;
  for (n = 3; n <= *nPS; ++n){
    ClearMat(pmat_next, *ndist, *ndist) ;
    ProbaNextConsistent(ISdist, pmat_next, pmat, *ndist, *IScn) ;
    probas[n-1] = MatSum(pmat_next, *ndist, *ndist) ;
    for (i = 0 ; i < *ndist; ++i)
      for (j = 0 ; j < *ndist; ++j)
	pmat[i][j] = pmat_next[i][j];
    //Rprintf("#") ;
  }
  //Rprintf("\n") ;
  Free(pmat); 
  Free(pmat_next) ;
}

void pdist_simu(int *dn, int *L, int *nPS, double * psim , int* nsim){
  int i, n;
  double cmin, cmax;
  double r ;
  double mu = ((double) *dn)/ (*L) ;
  int ncons = 0 ;
  GetRNGstate();
  for (i=0 ; i < *nsim; ++i){
    cmin = 1 ; 
    cmax = 0 ;
    for (n =0; n < *nPS; ++n){
      r = unif_rand();
      cmin = fmin2(r, cmin);
      cmax = fmax2(r, cmax);
    }
    //Rprintf("simu n %d: max :%.2f, min: %.2f\n", i, cmax , cmin) ;
    if ((cmax - cmin) < mu)
      ncons++  ;
  }
  PutRNGstate();
  *psim = ((double) ncons) / *nsim ;
}


/* void pIS_simu(double* ISdist, int *ndist, double *pest, int *nPS, int *IScn, int *nsim){ */
/*   int i, n;  */
/*   int c = 0 ; */
/*   int cmin, cmax; */
/*   time_t tt ; */
/*   //tt = time(NULL) ; */
/*   //set_seed(tt); */
/*   GetRNGstate(); */
/*   for (i = 0 ; i < *nsim; ++i){ */
/*     cmin = cmax = ISdist[(int) unif_rand() * (*ndist)]; */
/*     for (n = 1; n < *nPS; *nPS++) */
/*       { */
/* 	cmin = imin */
/*       } */
/*     if () */
/*       c++ */
/*   } */
/*   pest = ((double) c) / nsim ; */
/*   PutRNGstate(); */
/* }  */

