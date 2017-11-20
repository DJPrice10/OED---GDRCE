/* Function to evaluate the metric, sum(abs(x-y))/std(x) */
/* Written by David Price, 27/11/2014 */

/* Input alldata, observeddata, standard deviations */


#include "mex.h"
#include <math.h>



void discrepfcn(double all[], double obs[], double sd[], double disc[], int nr, int nc)
{
    int row, col;
    double discs[nr][nc];
    
    /* Evaluate discrepancies between each row and column 
     * between obs dat and all data
     */
        for(col=0; col<nc; col++){
            for(row=0; row<nr; row++){
                discs[row][col]=fabs(all[row+nr*col]-obs[col])/sd[col];
                disc[row]= disc[row]+discs[row][col]; 
            }
        }
    /* Sum each row 
        for(row=0; row<nr; row++)
            for(col=0; col<nc; col++)
                disc[row]= disc[row]+discs[row][col]; 
    */
}

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* double *alldat, *obsdat, *stds, *disc; */
  size_t mrows,ncols;
  /* size_t is an unsigned var type for dimensions */
  
  mrows = mxGetM(prhs[0]);
  ncols = mxGetN(prhs[0]);
  
  /* Create matrix for the return argument. */
  plhs[0] = mxCreateDoubleMatrix((mwSize)mrows, 1, mxREAL);
  
  /* Pass all data, observed data, std devs, and the pointer for result, plus row and col dims for defining discrepancy matrix.
   * Pass pointer to output to fcn (disc) so that disc is returned 
   */
  discrepfcn(mxGetPr(prhs[0]), mxGetPr(prhs[1]), mxGetPr(prhs[2]), mxGetPr(plhs[0]), mrows, ncols);

}
