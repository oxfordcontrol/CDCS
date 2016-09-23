/***********************************************************************
 * mexsmat.c : C mex file
 *
 *  M = smat(x);
 *
 *   Input: x     = n(n+1)/2 vector
 *
 * Based on code from SDPT3.4
 ***********************************************************************/

#include <mex.h>
#include <math.h>
#include <matrix.h>


/**********************************************************
 * form Q using the lower triangular part
 **********************************************************/
void sym(double *Q, int n)
{
    int j, k, jn, kn;
    
    for (j=0; j<n; ++j){
        jn = j*n;
        for (k=0; k<j ; ++k){
            kn = k*n;
            Q[k+jn] = Q[j+kn]; }
    }
    return;
}

/**********************************************************
 * form lower triangular part of P
 * single block
 **********************************************************/
void smat1(int n, const double ir2,
        const double *A, const mwIndex *irA, const mwIndex *jcA, int isspA,
        double *B, mwIndex *irB, mwIndex *jcB, int isspB)
        
{  int idx, i, j, r, jn, k, kstart, kend, idxj, j2, count;
   double tmp;
   double hf=0.5;
   
   /*dense input, dense output*/
   if (!isspA & !isspB) {
       idx = 0;
       for (j=0; j<n; j++) {
           jn = j*n;
           B[j+jn] = A[idx];
           idx++;
           for (i=j+1; i<n; i++) {
               B[i+jn] = ir2*A[idx];
               idx++; }
       }
   }
   
   /*sparse input, sparse output*/
   else if (isspA & isspB) {
       count = 0;
       j2 = 0; 
       idxj = 0;
       kstart = jcA[0];  
       kend = jcA[1];
       for (k=kstart; k<kend; k++) { /*for each nonzero entry of A*/
           r = irA[k];               /*which entry?*/
           
           /*Compute subscripts for linear index r*/
           for (j=j2; j<n; j++) {   
               i=r-idxj;             
               if (i>n-j-1) {            /*if i not a valid row index*/
                   idxj+= n-j;       /*update index*/
               }
               else {
                   break;           
               }
           } 
            /* which column? */
           j2=j;
           
           /*Set entry*/
           if (i==0) {
               irB[count] = i+j;
               B[count] = hf*A[k];
               ++jcB[j+1];
           }  
           else {
               irB[count] = i+j;
               B[count] = ir2*A[k];
               ++jcB[j+1];
           }
           count++;
       }
       
       /*cumsum of jcB*/
       for (j=0; j<n; j++) { 
           jcB[j+1] += jcB[j]; 
       }
   }
   
   if (!isspB) {
       sym(B,n);
   }
   return;
}

/**********************************************************
 * MAIN MEX FUNCTION
 ***********************************************************/
void mexFunction(int nlhs,   mxArray  *plhs[],
        int nrhs,   const mxArray  *prhs[] )
        
{    mxArray  *rhs[2];
     mxArray  *blk_cell_pr, *A_cell_pr;
     double   *A,  *B,  *blksize, *AI, *BI;
     mwIndex  *irA, *jcA, *irB, *jcB;
     int       mblk, mA, nA, m1, n1, rowidx, colidx, isspA, isspB;
     int      n, n2, k, nsub, index, numblk, NZmax;
     double   ir2=1/sqrt(2);
     
     /* CHECK FOR PROPER NUMBER OF ARGUMENTS */
     if (nrhs != 1){
         mexErrMsgTxt("smat: requires 1 input argument."); }
     else if (nlhs>1){
         mexErrMsgTxt("smat: requires 1 output argument."); }
     
     
     /***** assign pointers *****/
     A  = mxGetPr(prhs[0]);
     m1 = mxGetM(prhs[0]);
     n1 = mxGetN(prhs[0]);
     if (n1==1) {
         n = ( sqrt(8*m1+1)-1 )/2;
     }
     else if (m1==1) {
         n = ( sqrt(8*n1+1)-1 )/2;
     }
     else {
         mexErrMsgTxt("smat: is only defined for vector arguments!"); 
     }
     
     
     isspA = mxIsSparse(prhs[0]);
     if (isspA) {
         irA = mxGetIr(prhs[0]);
         jcA = mxGetJc(prhs[0]);
     }
     
     /***** create return argument *****/
     if (isspA) {
         NZmax = jcA[1]-jcA[0];
         rhs[0] = mxCreateSparse(n,n,2*NZmax,mxREAL);
         B = mxGetPr(rhs[0]);
         irB = mxGetIr(rhs[0]);
         jcB = mxGetJc(rhs[0]);
     }
     else {
         plhs[0] = mxCreateDoubleMatrix(n,n,mxREAL);
         B = mxGetPr(plhs[0]);
     }
     
     /***** Do the computations in a subroutine *****/
     smat1(n,ir2,A,irA,jcA,isspA,B,irB,jcB,isspA);
     if (isspA) {
         /*** if isspB, (actual B) = B+B' ****/
         mexCallMATLAB(1, &rhs[1], 1, &rhs[0], "ctranspose");
         mexCallMATLAB(1, &plhs[0],2, rhs, "+");
         mxDestroyArray(*rhs);
     }
     
     
     return;
}
/**********************************************************/

