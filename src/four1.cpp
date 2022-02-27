#include <stdio.h>
#include <math.h>
#include <fftw3.h>
#include <iostream>
#include "Utilities.hpp"

#define PI 3.141592653589793238462643383279

extern fftw_plan  p1_FFTW;
extern fftw_plan  p2_FFTW;


void four1_FFTW(Vektor<double>& data, unsigned long nn, int isign)
{

/*******************************************************************
*
*  This is a wrapper for the FFT function four1 in Numerical Recipes
*      It follows all of the argument specs there.
*
*  NB:
*      1. size of data is 2*nn, not nn
*      2. range of data is 1 to 2*nn. NOT 0 to 2*nn-1
*      3. if calling with data[0..2*nn-1], use data-1 in input argument
*         i.e. four1_FFTW(data-1, nn, isign);
*
*
* ** add these to the right places
*
*     #include <fftw3.h>               ** header
*
*
*     fftw_plan p1_FFTW;               **  global variables to store FFTW wisdom
*     fftw_plan p2_FFTW;
*
*
*     extern fftw_plan     p1_FFTW;    ** header of function file
*     extern fftw_plan     p2_FFTW;
*
*
*     fftw_destroy_plan(p1_FFTW);      ** end of main program
*     fftw_destroy_plan(p2_FFTW);
*
*
*  **  compile and link with
*
*     gcc -o main main.c four1_FFTW.c -L./ -llibfftw3-3 -lm
*
*
*  **  library files to add to directory
*
*       libfftw3-3.dll,  fftw3.h
*
*
********************************************************************/
  int n=nn, i;

  /* ****************************************************
   * Allocate memory
   * ****************************************************/
  fftw_complex *b1, *b2;

  b1 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);
  b2 =(fftw_complex*) fftw_malloc(sizeof(fftw_complex)*n);

  /* ****************************************************
   * Create forward FFT plan from b1 into b2
   * ****************************************************/
  p1_FFTW=fftw_plan_dft_1d(n, b1, b2, FFTW_FORWARD, FFTW_ESTIMATE);

  /* ****************************************************
   * Create reverse FFT plan from b1 into b2
   * ****************************************************/
  p2_FFTW=fftw_plan_dft_1d(n, b1, b2, FFTW_BACKWARD, FFTW_ESTIMATE);

  /* ****************************************************
   * FFT:  insert - sign for imag part because
   *       Num Rec FFT  uses exp(+i*2*pi*m*n/N) whereas
   *               FFTW uses exp(-i*2*pi*m*n/N)
   * ****************************************************/
    for (i=0; i < n; i++){
      b1[i][0] =  data(2*i+1);
      b1[i][1] = -data(2*i+2);
    }

    if (isign ==  1)  fftw_execute(p1_FFTW);    /* FFT             */
    if (isign == -1)  fftw_execute(p2_FFTW);    /* IFFT scale by n */

    for (i=0; i < n; i++){
      data(2*i+1) =  b2[i][0];
      data(2*i+2) = -b2[i][1];
    }

  /* **** Free memory   **************************/
  fftw_free(b1);
  fftw_free(b2);

}

