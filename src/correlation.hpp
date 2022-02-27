#ifndef CORRELATION_HPP_INCLUDED
#define CORRELATION_HPP_INCLUDED
#include "Utilities.hpp"
//#include "lineshape.hpp"
void correlation(Vektor<double> parameter_vec, Vektor<double>& CIm);
void PM_0_Four(int M, Vektor<double>& DMReals,Vektor<double> GReal, Vektor<double> GImag, double te, double ntime, Matrix<double> evecs, Vektor<double> evals, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat);

#endif // CORRELATION_HPP_INCLUDED
