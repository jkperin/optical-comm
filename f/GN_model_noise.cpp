/* Compute nonlinear noise power at each channel 
 * This is a C++ implementation of the Matalb function GN_model_noise_m.m. 
 * The C++ implementation runs about 100x faster than the Matlab 
 * implementation even when the Matlab implementation uses parallelism
 * 
 * Usage: 
 * > NL = GN_model_noise(P, D)
 * 
 * Inputs:
 * P is a 1 x N vector containing the Power at each of the N channels
 * D is a cell array of length 3, such that D{l+2} corresponds to the nonlinear 
 * coefficients for l = -1, 0, or 1. Each entry of D is a 2*N-1 x 2*N-1 matrix
 *
 * Output: 
 * NL is a 1 x N vector containing the nonlinear noise power at each channel.
 * NL is only computed for non-zero power channels
 */

#include "mex.h"
#include "matrix.h"

void GN_model_noise(double * NL, double * P, double * Dn, double * D0, double * Dp, int N, int M)
{   
    int idx = 0, i = 0, j = 0, ij = 0;
    int Ncenter = (M+1)/2; // index to the center of D matrices
    for(int n = 0; n < N; n++) {
        if (P[n] == 0) 
            continue;
        NL[n] = 0;
        for(int n1 = 0; n1 < N; n1++) {
            if (P[n1] == 0) 
                continue;
            for(int n2 = 0; n2 < N; n2++) {
                if (P[n2] == 0) 
                    continue;

                i = Ncenter - (n1 - n) - 1;
                j = Ncenter + (n2 - n) - 1;
                ij = j*M + i; // indexing of multdimensional vector as a 1 x M dimensional vector 
                // Note: MATLAB 2D matrix memory is arranged in column-order
                // (like Fortran), but C 2D array memory is arranged in row-order. 
                
                /* Compute nonlinear noise power for different values of l */
                // l = -1
                idx =  n1 + n2 - n - 1 - 1;
                if (idx >= 0 && idx < N && P[idx] != 0)
                    NL[n] += P[n1]*P[n2]*P[idx]*Dn[ij];
                
                // l = 0
                idx =  n1 + n2 - n - 1;
                if (idx >= 0 && idx < N && P[idx] != 0)
                    NL[n] += P[n1]*P[n2]*P[idx]*D0[ij];
                
                // l = 1
                idx =  n1 + n2 - n -1 + 1;
                if (idx >= 0 && idx < N && P[idx] != 0)
                    NL[n] += P[n1]*P[n2]*P[idx]*Dp[ij];       
            }
        }
    }
    
    return;
}

/*************************************************************************/
/************************* GATEWAY ROUTINE *******************************/
/*************************************************************************/

// NL = GN_model_noise(P, D)
// P is a vector of length N
// D is a cell array of length 3. Each entry of D is a square matrix of size 2*N-1 x 2N-1.

// nlhs	Number of output (left-side) arguments, or the size of the plhs array.
// plhs	Array of output arguments.
// nrhs	Number of input (right-side) arguments, or the size of the prhs array.
// prhs	Array of input arguments.

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
    mxArray *mxDn, *mxD0, *mxDp;
    double *P, *Dn, *D0, *Dp; // inputs
    double * NL; // output
    int N; // number of channels N = length(P)
    int M; // number of rows in nolinear coefficient matrix. M > 2*N-1
      
    if(nrhs != 2)
        mexErrMsgTxt("GN_model_noise takes exactly two arguments.");
    if(nlhs != 1)
        mexErrMsgTxt("GN_model_noise has only one output");
    
    P = mxGetPr(prhs[0]);
    N = mxGetN(prhs[0]);
    
    // check if inputs are consistent
    if (mxGetM(prhs[0]) != 1 || !mxIsDouble(prhs[0]))
        mexErrMsgTxt("First argument must be a vector.");
   
    size_t dims = mxGetNumberOfElements(prhs[1]);
    if (!mxIsCell(prhs[1]) || dims != 3)
        mexErrMsgTxt("Second argument must be a cell array of length 3.");
     
    mxDn = mxGetCell(prhs[1], 0);
    mxD0 = mxGetCell(prhs[1], 1);
    mxDp = mxGetCell(prhs[1], 2);
    
    if (mxGetN(mxDn) != mxGetM(mxDn) || mxGetN(mxD0) != mxGetM(mxD0) || mxGetN(mxDp) != mxGetM(mxDp))
        mexErrMsgTxt("Elements of cell array must be square matrices");
    
    if (mxGetN(mxDn) < 2*N-1 || mxGetN(mxD0) < 2*N-1 || mxGetN(mxDp) < 2*N-1)
        mexErrMsgTxt("Matrices in cell array must be at least 2*N-1 x 2*N -1, where N is the length of the vector passed as first argument");
    
    M = mxGetN(mxDn);
    Dn = mxGetPr(mxDn);
    D0 = mxGetPr(mxD0);
    Dp = mxGetPr(mxDp);
    
    plhs[0] = mxCreateDoubleMatrix(1, N, mxREAL);
    NL = mxGetPr(plhs[0]);
    
    // GN_model_noise(double * NL, double * P, double * Dn, double * D0, double * Dp, int N, int M)
    GN_model_noise(NL, P, Dn, D0, Dp, N, M);
    
    return;
}