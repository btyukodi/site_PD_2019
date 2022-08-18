#ifndef FFT2D_H
#define FFT2D_H
/********************************************************************
   F A S T   F O U R I E R   T R A N S F O R M   P R O G R A M S

  by Wang Jian-Sheng 4 Nov 1998, added fft2D(), 11 Apr 2003 
---------------------------------------------------------------------

  Reference: "Computational Frameworks for the Fast Fourier 
              Transform", Charles Van Loan, SIAM, 1992.

  There are many FFT algorithms, the most important ones are
     COOLEY-TUKEY:  in place, bit reversal
     STOCKHAM AUTOSORT:  additional memory size of input data
     MIXED RADIX:  20% less operations comparing to Cooley-Tukey
     PRIME FACTOR: arbitrary length n

  We use a combination of the Stockham autosort algorithm 1.7.2, 
  page 57, and multirow Cooley-Tukey (3.1.7), page 124, of the 
  reference above.  

  The discrete Fourier transform is defined by
  y[k] = sum_(j=0,n-1) x[j] exp(-2 Pi sqrt[-1] j k/n), 
  k=0,1,...,n-1.  The factor (1/n) is not included.  
  If y[]<-x[]; fft(x,n,1); fft(x,n,-1); then y[]==x[]/n is true.
  Three dimensional transform is generalized straightforwardly.

   Interface and usage:
   1D Fourier transform 
   Use: fft(x, n, flag)
      x    : an array of structure type complex;
      n    : the size of data, must be a power of 2;
      flag : 1 for forward transform, -1 for inverse transform.

   3D Fourier transform
   Use :  fft3D(x, n1, n2, n3, flag)
     x    : 1D array of type complex representing 3D array; 
            mapping through C convention, i.e., 
            (i,j,k) -> k + n3*j + n2*n3*i;
     n1, n2, n3 : dimensions in three directions;
     flag : same as in 1D.

   2D FFT is similar but with n1 and n2 only.

**********************************************************************/
  
/* Data type and new names for flexibility:
  
    real:    Basic data type for floating point computations
             (typedef double  real;)
    complex: Structure for complex numbers, real and imaginary parts 
             are referred as c.Re, c.Im.
             (typedef struct { real Re; real Im; }  complex;)
*/
/* Inclusion of standard C libraries */


typedef double real;                 /* can be long double, double, or float */
typedef struct { real Re; real Im; }  complex; 

void fft(complex x[], int n, int flag);
void fft2D(complex x[], int n1, int n2, int flag);
void fft3D(complex x[], int n1, int n2, int n3, int flag);

#endif
