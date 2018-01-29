/* Compute gradient of nonlinear noise power at each channel 
 * Usage: 
 * > [dNL, NL] = GN_model_noise_gradient(P, D)
 * 
 * Inputs:
 * P is a 1 x N vector containing the Power at each of the N channels
 * D is a cell array of length 3, such that D{l+2} corresponds to the nonlinear 
 * coefficients for l = -1, 0, or 1. Each entry of D is a 2*N-1 x 2*N-1 matrix
 *
 * Output: 
 * dNL is a 1 x N vector containing the nonlinear noise power gradient at each channel.
 * dNL is indexed such that
 *          dNL(i, j) = d NL_j / d P_i
 * NL is a 1 x N vector containing the nonlinear noise power at each channel.
 * NL is only computed for non-zero power channels
 */

#include "mex.h"
#include "matrix.h"

void GN_model_noise_gradient(double * dNL, double * NL, double * P, double * Dn, double * D0, double * Dp, int N)
{   
    int idx = 0, i = 0, j = 0, ij = 0;
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

                i = N - (n1 - n) - 1;
                j = N + (n2 - n) - 1;
                ij = j*(2*N-1) + i; // indexing of multdimensional vector as a 1 x N dimensional vector 
                // Note: MATLAB 2D matrix memory is arranged in column-order
                // (like Fortran), but C 2D array memory is arranged in row-order. 
                
                /* Compute nonlinear noise power for different values of l */
                // l = -1
                idx =  n1 + n2 - n - 1 - 1;
                if (idx >= 0 && idx < N && P[idx] != 0) {
                    NL[n] += P[n1]*P[n2]*P[idx]*Dn[ij];
                                
                    for (int k = 0; k < N; k++) {
                        if (k == n1 && k != n2 && k != idx)
                            dNL[n*N + k] += P[n2]*P[idx]*Dn[ij];
                        else if (k != n1 && k == n2 && k != idx)
                            dNL[n*N + k] += P[n1]*P[idx]*Dn[ij];
                        else if (k != n1 && k != n2 && k == idx)
                            dNL[n*N + k] += P[n1]*P[n2]*Dn[ij];
                        else if (k == n1 && k == n2 && k != idx) // 2 equal terms
                            dNL[n*N + k] += 2*P[n1]*P[idx]*Dn[ij];
                        else if ((k != n1 && k == n2 && k == idx) || (k == n1 && k != n2 && k == idx))
                            dNL[n*N + k] += 2*P[n2]*P[n1]*Dn[ij];
                        else if (k == n1 && k == n2 && k == idx) // all terms are equal
                            dNL[n*N + k] += 3*P[n1]*P[n2]*Dn[ij];
                    }
                }

                // l = 0
                idx =  n1 + n2 - n - 1;
                if (idx >= 0 && idx < N && P[idx] != 0) {
                    NL[n] += P[n1]*P[n2]*P[idx]*D0[ij];
                    
                    for (int k = 0; k < N; k++) {
                        if (k == n1 && k != n2 && k != idx)
                            dNL[n*N + k] += P[n2]*P[idx]*D0[ij];
                        else if (k != n1 && k == n2 && k != idx)
                            dNL[n*N + k] += P[n1]*P[idx]*D0[ij];
                        else if (k != n1 && k != n2 && k == idx)
                            dNL[n*N + k] += P[n1]*P[n2]*D0[ij];
                        else if (k == n1 && k == n2 && k != idx) // 2 equal terms
                            dNL[n*N + k] += 2*P[n1]*P[idx]*D0[ij];
                        else if ((k != n1 && k == n2 && k == idx) || (k == n1 && k != n2 && k == idx))
                            dNL[n*N + k] += 2*P[n2]*P[n1]*D0[ij];
                        else if (k == n1 && k == n2 && k == idx) // all terms are equal
                            dNL[n*N + k] += 3*P[n1]*P[n2]*D0[ij];
                    }
                }
                
                // l = 1
                idx =  n1 + n2 - n -1 + 1;
                if (idx >= 0 && idx < N && P[idx] != 0) {
                    NL[n] += P[n1]*P[n2]*P[idx]*Dp[ij];  
                    
                    for (int k = 0; k < N; k++) {
                        if (k == n1 && k != n2 && k != idx)
                            dNL[n*N + k] += P[n2]*P[idx]*Dp[ij];
                        else if (k != n1 && k == n2 && k != idx)
                            dNL[n*N + k] += P[n1]*P[idx]*Dp[ij];
                        else if (k != n1 && k != n2 && k == idx)
                            dNL[n*N + k] += P[n1]*P[n2]*Dp[ij];
                        else if (k == n1 && k == n2 && k != idx) // 2 equal terms
                            dNL[n*N + k] += 2*P[n1]*P[idx]*Dp[ij];
                        else if ((k != n1 && k == n2 && k == idx) || (k == n1 && k != n2 && k == idx))
                            dNL[n*N + k] += 2*P[n2]*P[n1]*Dp[ij];
                        else if (k == n1 && k == n2 && k == idx) // all terms are equal
                            dNL[n*N + k] += 3*P[n1]*P[n2]*Dp[ij];
                    }
                }
            }
        }
    }
    
    return;
}

/*************************************************************************/
/************************* GATEWAY ROUTINE *******************************/
/*************************************************************************/

// [dNL, NL] = GN_model_noise_gradient(P, D)
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
    double * dNL, * NL; // output
    int N; // number of channels N = length(P)
      
    if(nrhs != 2)
        mexErrMsgTxt("GN_model_noise takes exactly two arguments.");
    if(nlhs != 2)
        mexErrMsgTxt("GN_model_noise takes exactly two outputs");
    
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
    
    if (mxGetN(mxDn) != 2*N-1 || mxGetN(mxD0) != 2*N-1 || mxGetN(mxDp) != 2*N-1)
        mexErrMsgTxt("P and D dimensions are not consistent");
    
    Dn = mxGetPr(mxDn);
    D0 = mxGetPr(mxD0);
    Dp = mxGetPr(mxDp);
    
    // Create outputs
    plhs[0] = mxCreateDoubleMatrix(N, N, mxREAL);
    dNL = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(1, N, mxREAL);
    NL = mxGetPr(plhs[1]);
    
    // GN_model_noise_gradient(double * dNL, double * NL, double * P, double * Dn, double * D0, double * Dp, int N)
    GN_model_noise_gradient(dNL, NL, P, Dn, D0, Dp, N);
    
    return;
}