/***********************************************************************
 * mexsvec.c : C mex file
 *
 *   x = svec(A);
 *
 *   Input: A     = nxn matrix
 *
 * Based on code from SDPT3.4
 ***********************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>


/**********************************************************
 * Stack lower triangular part of A row-wise into a column vector
 **********************************************************/
void svec1(int n, double r2,
        double *A, mwIndex *irA, mwIndex *jcA, int isspA,
        double *B, mwIndex *irB, mwIndex *jcB, int isspB)
        
{  int idx, rowidx, i, j, jn, k, kstart, kend;
   
   /* Dense input, dense output */
   if (!isspB & !isspA) {
       idx = 0;
       for (j=0; j<n; j++) {
           B[idx] = A[j+j*n];
           idx++;
           for (i=j+1; i<n; i++) {
               B[idx] = A[i+n*j]*r2;
               idx++;
           }
       }
       
   }
   
   /* Sparse input, sparse output */
   else if (isspB & isspA) {
       idx = 0;
       for (j=0; j<n; j++) {    /*for each column*/
           kstart = jcA[j];     /*how many nonzero els so far*/
           kend = jcA[j+1];     /*how many nonzero els after col j is counted*/
           for (k=kstart; k<kend; k++) { /*for each nonzero entry of col j*/
               i = irA[k];      /*the row index of kth nnz element of col j*/
               if (i > j) {
                   irB[idx] = i + n*j - ((j+1)*j)/2;
                   B[idx] = A[k]*r2;
                   idx++; }
               else if (i==j) {
                   irB[idx] = i + n*j - ((j+1)*j)/2;
                   B[idx] = A[k];
                   idx++;
               }
           }
       }
       jcB[1] = idx;
   }
   return;
}

/**********************************************************
 * Stack lower triangular part of A row-wise into a column vector
 **********************************************************/
void copyvec(int n,
        double *A, mwIndex *irA, mwIndex *jcA, int isspA,
        double *B, mwIndex *irB, mwIndex *jcB, int isspB)
        
{  int idx, rowidx, i, j, jn, k, kstart, kend;
   
   /* Dense input, dense output */
   if (!isspB & !isspA) {
       idx = 0;
       for (j=0; j<n; j++) {
           B[idx] = A[j];
           idx++;
       }
   }
   /* Sparse input, sparse output */
   else if (isspB & isspA) {
       idx = 0;
       kstart = jcA[0];   /*how many nonzero els so far*/
       kend = jcA[1];     /*how many nonzero els after col j is counted*/
       for (k=kstart; k<kend; k++) { /*for each nonzero entry of col j*/
           irB[idx] = irA[k];      /*the row index of kth nnz element of col j*/
           B[idx] = A[k];
           idx++;
       }
       jcB[1] = idx;
   }
   
   return;
   
}


/**********************************************************
 * MAIN MEX FUNCTION
 ***********************************************************/
void mexFunction(int nlhs, mxArray  *plhs[],
        int nrhs, const mxArray  *prhs[] )
        
{    double   *A,  *B,  *AI, *BI, *blksize;
     mwIndex  *irA, *jcA, *irB, *jcB;
     int      isspA, isvecA;
     int      m, n, n2, k, index, NZmax;
     double   r2 = sqrt(2);;
     
     /* CHECK FOR PROPER NUMBER OF ARGUMENTS */
     
     if (nrhs != 1){
         mexErrMsgTxt("svec: requires 1 input argument."); }
     if (nlhs > 1){
         mexErrMsgTxt("svec: requires 1 output argument."); }
     
     /* CHECK THE DIMENSIONS */
     
     m = mxGetM(prhs[0]);
     n = mxGetN(prhs[0]);
     if (n==1) {
         isvecA = 1;
         n2 = m;
     }
     else if (m == n) {
         isvecA = 0;
         /* Find size of svectorized variable */
         n2 = n*(n+1)/2;
     }
     else {
         mexErrMsgTxt("svec: matrix must be square.");
     }
     
     
     
     
     /***** assign pointers *****/
     A = mxGetPr(prhs[0]);
     isspA = mxIsSparse(prhs[0]);
     if (isspA) {
         irA = mxGetIr(prhs[0]);       /* get row indices */
         jcA = mxGetJc(prhs[0]);       /* get col indices */
         NZmax = mxGetNzmax(prhs[0]);  /* # allocated nonzero elements */
     }
     else {
         NZmax = n2;
     }
     
     
     /***** create return argument *****/
     if (isspA) {
         plhs[0] = mxCreateSparse(n2,1,n2,mxREAL);
         B = mxGetPr(plhs[0]);
         irB = mxGetIr(plhs[0]);
         jcB = mxGetJc(plhs[0]);
         jcB[0] = 0;
     }
     else {
         plhs[0] = mxCreateDoubleMatrix(n2,1,mxREAL);
         B = mxGetPr(plhs[0]);
     }
     
     /***** Do the computations in a subroutine *****/
     if (isvecA) {
         copyvec(m,A,irA,jcA,isspA,B,irB,jcB,isspA);
     }
     else {
         
         svec1(n,r2,A,irA,jcA,isspA,B,irB,jcB,isspA);
     }
     return;
}
/**********************************************************/
