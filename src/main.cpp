#include <bits/stdc++.h>
#include "init.hpp"
#include "excmatrix.hpp"
#include "randomgauss.hpp"
#include "correlation.hpp"
#include "Utilities.hpp"
#include "excdomains.hpp"
#include <fftw3.h>
#include "LU_Decomp.hpp"

using namespace std;

int main(){

    Vektor<double> debye_vec;
    Vektor<double> sitevec;
    Vektor<double> parameter_vec;
    Vektor<string> chlname_vec;
    Vektor< Vektor< Vektor<double> > > dipolemat;
    Matrix<double> evecs;
    Vektor<double> evals;

    Vektor<double> Huang_Fac;
    Vektor<double> Inh_Broads;
    Matrix<double> Exc2Shifts;

    //-----Get init.dat file and dipole vectors from nb_nb.dat-----
	init(debye_vec, parameter_vec, Huang_Fac, Inh_Broads);
	getdipoles(dipolemat, debye_vec, chlname_vec);

	//-----Get site-energies from Site_Energies.dat file as well as the excitonic couplings from excitonic_couplings.dat-----
    getSiteEnergies(sitevec, chlname_vec);
    diagonalizeexcmatrix(evecs, evals, chlname_vec, parameter_vec, sitevec);
    get2Exc_Shifts(Exc2Shifts, sitevec);

    //-----Setup the Helix-----
    Vektor< ChlAtoms > Chl_Vec;
//    setup_Helix(Chl_Vec, dipolemat);
    setup_Helix_2StateModel(Chl_Vec, dipolemat);

//    Lin_Pump_Test(sitevec, parameter_vec, chlname_vec, dipolemat);
    CD_Pump_Test(sitevec, parameter_vec, chlname_vec, dipolemat, Exc2Shifts);

    return 0;
}
