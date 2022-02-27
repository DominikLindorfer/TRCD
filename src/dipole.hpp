#ifndef DIPOLE_HPP_INCLUDED
#define DIPOLE_HPP_INCLUDED
#include "Utilities.hpp"
#include "ChlClass.hpp"

using namespace std;
void getdipoles(Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double> debye_vec, Vektor<string>& chlname_vec);
void get_nanc_nbnd(Vektor< Vektor< Vektor<double> > >& dipolemat);
void turn_dipoles(Vektor< Vektor< Vektor<double> > >& dipolemat, double t_ang);
void get_coordinates(Vektor< ChlAtoms >& Chl_Vec, string FileName);
int get_charges(ChlAtoms& Chl_Vec, string FileName);
int get_BChl_charges(ChlAtoms& Chl_Vec, string FileName);

#endif // DIPOLE_HPP_INCLUDED
