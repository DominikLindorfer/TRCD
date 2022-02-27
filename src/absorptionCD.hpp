#ifndef ABSORPTIONCD_HPP_INCLUDED
#define ABSORPTIONCD_HPP_INCLUDED

#include <fftw3.h>
#include "Utilities.hpp"
#include <iostream>
#include <fstream>
#include "correlation.hpp"
#include "ChlClass.hpp"
#include <vector>

#include "lineshape.hpp"

using namespace std;

void crossprod(Vektor<double> A,Vektor<double> B,Vektor<double>& R){
    //because uBLAS is for d-dimensional vectors
    R(0)= A(1)*B(2) - B(1)*A(2);
    R(1)= B(0)*A(2) - A(0)*B(2);
    R(2)= A(0)*B(1) - A(1)*B(0);
}

Vektor<double> Cross(Vektor<double> A, Vektor<double> B){
    Vektor<double> R(3);

	R(0)= A(1)*B(2) - B(1)*A(2);
    R(1)= B(0)*A(2) - A(0)*B(2);
    R(2)= A(0)*B(1) - A(1)*B(0);

    return R;
}

void CDabsorption_domains(Matrix<double> evecs, Vektor<double> evals, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor<double>& CD, Vektor<double>& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p){

    //Allocate a RM Vector for all RM values of eigenvalues M
    Vektor<double> RM(evals.size());
    for(int i=0;i < (int)RM.size();i++){
        RM(i)=0;
    }

    Vektor<double> CrossResult(3);

    //ABSORPTION-PART
    Vektor< Vektor<double> > trans_dipole_moment(evals.size());
    ublas::zero_vector<double> zero(3);
    for(int i=0;i < (int)trans_dipole_moment.size();i++){
        trans_dipole_moment(i)=zero;    //zerovectors so += can be used!
    }

    for(int i = 0; i < (int)trans_dipole_moment.size(); i++){              //check if error occurs!!
        for(int j = 0; j < (int)trans_dipole_moment.size(); j++){
            trans_dipole_moment(i)+=dipolemat(0)(j)*evecs(i,j);      //careful! dipolemat(0) because it contains the ei_dipole
        }                                                           //there are only 2 entries until now!!!
    }

    //cout<<"Dipolemat"<<dipolemat(0)<<endl;     //errorchecking
//    cout<<"Transdipolemoment: "<<trans_dipole_moment<<endl;
    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size() ;M++){
        for(int i=0;i < (int)evals.size()-1;i++){
            for(int j=i+1;j < (int)evals.size();j++){

                //crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);
                crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);

                //crossprod(dipolemat(0)(0),dipolemat(0)(1),CrossResult);
                //dipolemat(1,M) because entry 1 contains the Centers of CHL

                RM(M) += 1.0 * evecs(M,i) * evecs(M,j) * inner_prod( dipolemat(1)(i) - dipolemat(1)(j) , CrossResult);
//                cout<<"i: "<<i<<"j: "<<j<<endl;   //errorchecking
//                cout << "M = " << M << " i,j " << i << j << "dipolemat(i) = " << dipolemat(0)(i) << "  dipolemat(j) " << dipolemat(0)(j) << endl;
//                cout << "evecs(M,i) = " << evecs(M,i) << "  evecs(M,j)= " << evecs(M,j) << endl;
//                cout << "dipolemat(1)(i) = " << dipolemat(1)(i) << " dipolemat(1)(j) = " << dipolemat(1)(j) << endl;
//                cout << "RM (" << M <<") = " << RM(M) << endl;
            }
        }

        //CHANGE LINESHAPEFUNCTION DM_WZ HERE!
        DM_WZ(M,DMReals,GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,CIm,evecs,evals,parameter_vec,dipolemat, in, out); //original 8192.0/16.0 ,8192.0 * 512.0  //Fast Version: 8192.0/256.0 ,8192.0 * 32.0 //NormalVersion 23.06.2015 8192.0/512.0 ,8192.0 * 16.0

        CD.resize(DMReals.size());
        CD = CD + DMReals * RM(M);

        //ABSORPTIONPART
        alpha.resize(DMReals.size());
        alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
        //cout<<"M = "<<M<<"  mu-M = "<<inner_prod(trans_dipole_moment(M),trans_dipole_moment(M))<<endl;
    }
}

void CDabsorption_domains_Get_Contributions(Matrix<double> evecs, Vektor<double> evals, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor< Vektor<double> >& CD, Vektor< Vektor<double> >& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p){

    //Allocate a RM Vector for all RM values of eigenvalues M
    Vektor<double> RM(evals.size());
    for(int i=0;i < (int)RM.size();i++){
        RM(i)=0;
    }

    Vektor<double> CrossResult(3);
	int LE = (int)evals.size() - 1;
	int HE = 0;

    //ABSORPTION-PART
    Vektor< Vektor<double> > trans_dipole_moment(evals.size());
    ublas::zero_vector<double> zero(3);
    for(int i=0;i < (int)trans_dipole_moment.size();i++){
        trans_dipole_moment(i)=zero;    //zerovectors so += can be used!
    }

    for(int i = 0; i < (int)trans_dipole_moment.size(); i++){              //check if error occurs!!
        for(int j = 0; j < (int)trans_dipole_moment.size(); j++){
            trans_dipole_moment(i)+=dipolemat(0)(j)*evecs(i,j);      //careful! dipolemat(0) because it contains the ei_dipole
        }                                                           //there are only 2 entries until now!!!
    }

    //cout<<"Dipolemat"<<dipolemat(0)<<endl;     //errorchecking
    //cout<<"Transdipolemoment: "<<trans_dipole_moment<<endl;
    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size() ;M++){
        for(int i=0;i < (int)evals.size()-1;i++){
            for(int j=i+1;j < (int)evals.size();j++){

                //crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);
                crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);

                //crossprod(dipolemat(0)(0),dipolemat(0)(1),CrossResult);
                //dipolemat(1,M) because entry 1 contains the Centers of CHL

                RM(M) += 1.0 * evecs(M,i) * evecs(M,j) * inner_prod( dipolemat(1)(i) - dipolemat(1)(j) , CrossResult);
//                cout<<"i: "<<i<<"j: "<<j<<endl;   //errorchecking
//                cout << "M = " << M << " i,j " << i << j << "dipolemat(i) = " << dipolemat(0)(i) << "  dipolemat(j) " << dipolemat(0)(j) << endl;
//                cout << "evecs(M,i) = " << evecs(M,i) << "  evecs(M,j)= " << evecs(M,j) << endl;
//                cout << "dipolemat(1)(i) = " << dipolemat(1)(i) << " dipolemat(1)(j) = " << dipolemat(1)(j) << endl;
//                cout << "RM (" << M <<") = " << RM(M) << endl;
            }
        }

        //CHANGE LINESHAPEFUNCTION DM_WZ HERE!
        DM_WZ(M,DMReals,GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,CIm,evecs,evals,parameter_vec,dipolemat, in, out); //original 8192.0/16.0 ,8192.0 * 512.0  //Fast Version: 8192.0/256.0 ,8192.0 * 32.0 //NormalVersion 23.06.2015 8192.0/512.0 ,8192.0 * 16.0

        //-----CD Part-----
        CD(0).resize(DMReals.size());
        CD(0) = CD(0) + DMReals * RM(M);

        //-----ABSORPTION PART-----
        alpha(0).resize(DMReals.size());
        alpha(0) = alpha(0) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));

        if(M == LE){
			//-----Low Energy Contribution-----
			CD(1).resize(DMReals.size());
			CD(1) = CD(1) + DMReals * RM(M);

			alpha(1).resize(DMReals.size());
			alpha(1) = alpha(1) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
        }

        else if(M == HE){
			//-----High Energy Contribution-----
			CD(2).resize(DMReals.size());
			CD(2) = CD(2) + DMReals * RM(M);

			alpha(2).resize(DMReals.size());
			alpha(2) = alpha(2) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
        }

    }
}

void CD_OD_Markov(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& CD, Vektor<double>& alpha, Vektor<double>& CIm, Vektor<Vektor<double> >& trans_dipole_moment){

	//-----Calculate CD and OD Spectrum with Markov Lineshape-----

    Vektor<double> RM(evals.size());
    for(int i = 0; i < (int)RM.size(); i++){
        RM(i) = 0;
    }

    Vektor<double> CrossResult(3);
    Vektor<double> DMReals;

	for(int M = 0; M < (int) RM.size(); M++){
		for(int i = 0; i < (int) evals.size() - 1; i++){
			for(int j = i + 1; j < (int) evals.size(); j++){

                crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);
                RM(M) += 1.0 * evecs(M,i) * evecs(M,j) * inner_prod( dipolemat(1)(i) - dipolemat(1)(j) , CrossResult);
            }
        }

        //-----Lineshape Function-----
        //DM_WZ(M,DMReals,GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,CIm,evecs,evals,parameter_vec,dipolemat, in, out);
        DM_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);

        //-----Result-----
        CD.resize(DMReals.size());
        CD = CD + DMReals * RM(M);

        alpha.resize(DMReals.size());
        alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
    }
}

void get_RotStrength(Matrix<double> evecs, Vektor<double> evals, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor<double>& RotStrength_vec, double E0){

    //Allocate a RM Vector for all RM values of eigenvalues M
    Vektor<double> RM(evals.size());
    for(int i=0;i < (int)RM.size();i++){
        RM(i)=0;
    }

    Vektor<double> CrossResult(3);

    //Create the exciton transition dipole moments
    Vektor< Vektor<double> > trans_dipole_moment(evals.size());
    ublas::zero_vector<double> zero(3);
    for(int i=0; i < (int)trans_dipole_moment.size(); i++){
        trans_dipole_moment(i)=zero;    //zerovectors so += can be used!
    }

    for(int i = 0; i < (int)trans_dipole_moment.size(); i++){
        for(int j = 0; j < (int)trans_dipole_moment.size(); j++){
            trans_dipole_moment(i)+=dipolemat(0)(j)*evecs(i,j);      //careful! dipolemat(0) because it contains the ei_dipole
        }
    }

    //cout<<"Dipolemat"<<dipolemat(0)<<endl;     //errorchecking
    //cout<<"Transdipolemoment: "<<trans_dipole_moment<<endl;
    Vektor<double> DMReals;

    for(int M = 0; M <  (int)RM.size(); M++){
        for(int i=0; i < (int)evals.size()-1; i++){
            for(int j=i+1; j < (int)evals.size(); j++){

                crossprod(dipolemat(0)(i),dipolemat(0)(j),CrossResult);
                //-----E0 [cm^-1] ; Rab [nm] ; mu [D] -----
                RM(M) += 2 * 1.7 * 1e-5 * E0 * evecs(M,i) * evecs(M,j) * 1e-1 * inner_prod( dipolemat(1)(i) - dipolemat(1)(j) , CrossResult) ;
//                cout << "M = " << M << " i,j " << i << j << "dipolemat(i) = " << dipolemat(0)(i) << "  dipolemat(j) " << dipolemat(0)(j) << endl;
//                cout << "evecs(M,i) = " << evecs(M,i) << "  evecs(M,j)= " << evecs(M,j) << endl;
//                cout << "dipolemat(1)(i) = " << dipolemat(1)(i) << " dipolemat(1)(j) = " << dipolemat(1)(j) << endl;
//                cout << "RM (" << M <<") = " << RM(M) << endl;
//                cout << "evals: " << evals(i) << "  " << evals(j) << endl;

            }
        }
        //cout<<"M = "<<M<<"  mu-M = "<<inner_prod(trans_dipole_moment(M),trans_dipole_moment(M))<<endl;
    }

    RotStrength_vec.resize(2);
    RotStrength_vec = RotStrength_vec + RM;

}

double dis_RotStrength(Vektor< Vektor<double> > RotStrength_sum, int side){
	//-----Calculate the <R->_dis (negative part)
	double sum = 0;

	for(int i = 0; i < (int)RotStrength_sum.size(); i++){

		if(RotStrength_sum(i)(0) > 0){
			cout << "Error in dis_RotStrength()" << endl;
		}

		sum += RotStrength_sum(i)(1);
	}

	return sum / (int)RotStrength_sum.size();
}

double W_KM(int K, int M, double delay, Vektor< Vektor<double> >& ci_evecs, Vektor<double>& evals_kinetic, Vektor<double>& P_eq, Vektor<double>& P0){
	//-----Helper-function W_KM(tau) which is related to PM(tau)-----

	double W_KM_t = 0;

	for(int i = 0; i < (int)evals_kinetic.size(); i++){
		W_KM_t += ci_evecs(i)(K) * sqrt(P_eq(K)) * ci_evecs(i)(M) * sqrt(P_eq(M)) * exp(delay * evals_kinetic(i));
	}

	return W_KM_t * P0(K) / P_eq(K);
}

double W_KM(int K, int M, double delay, Vektor< Vektor<double> >& ci_evecs, Vektor<double>& evals_kinetic, Vektor<double>& P_eq, Vektor<double>& P0, Vektor<Vektor<double> >& trans_dipole_moment){
	//-----Helper-function W_KM(tau) but normalized to |mu_K|-----

	double W_KM_t = 0;

	for(int i = 0; i < (int)evals_kinetic.size(); i++){
		W_KM_t += ci_evecs(i)(K) * sqrt(P_eq(K)) * ci_evecs(i)(M) * sqrt(P_eq(M)) * exp(delay * evals_kinetic(i));
	}

	return W_KM_t * P0(K) / P_eq(K) / inner_prod(trans_dipole_moment(K), trans_dipole_moment(K));
}

void Lin_GSBleach(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& alpha, Vektor<double>& CIm, Vektor<double>& GReals, Vektor<double>& GImags, fftw_complex* in, fftw_complex* out, Vektor<double>& P0, Vektor<Vektor<double> >& trans_dipole_moment){

    //-----Groundstate Bleaching-----

	Vektor<double> DMReals;

	for (int M = 0; M < (int)evals.size(); M++) {

		double P_0_sum = 0.0;

		//-----Sum of Exciton Populations at td = 0-----
		for(int K = 0; K < (int)evals.size(); K++){

			P_0_sum += P0(K);
		}

		//-----Lineshape Function-----
//		DM_WZ(M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);
		DM_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);

		//-----Result-----
		alpha.resize(DMReals.size());
		alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * P_0_sum * -1.0 / 9.0;
	}
}

void Lin_StimEmission(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& alpha, Vektor<double>& CIm, Vektor<double>& GReals, Vektor<double>& GImags, fftw_complex* in, fftw_complex* out, Vektor<double> P0, Vektor<Vektor<double> >& trans_dipole_moment, Vektor<double>& PM_delay){

    //-----Stimulated Emission-----

	Vektor<double> DMReals;

	//-----PURE DEPHASING-----
	//double puredeph = 1.0/0.1885;

	for (int M = 0; M < (int)evals.size(); M++) {

		//-----Lineshape Function-----
//		DM_Emission_WZ(M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);
		DM_Emission_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);
		//-----Result-----
		alpha.resize(DMReals.size());
		alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * PM_delay(M) * -1.0 / 9.0;
	}
}

void Lin_StimEmission_WKM(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& alpha, Vektor<double>& CIm, Vektor<double>& GReals, Vektor<double>& GImags, fftw_complex* in, fftw_complex* out, Vektor<Vektor<double> >& trans_dipole_moment, Vektor<double>& P0, Vektor<double>& P_eq, Vektor<double>& evals_kinetic, Vektor< Vektor<double> >& ci_evecs, double delay){

    //-----Stimulated Emission-----

	Vektor<double> DMReals;

	//-----PURE DEPHASING-----
	//double puredeph = 1.0/0.1885;

	for (int M = 0; M < (int)evals.size(); M++){
		//-----Lineshape Function-----
		DM_Emission_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);

		double W_KM_sum = 0;
		for(int K = 0; K < (int)evals.size(); K++){
			W_KM_sum += W_KM(K, M, delay, ci_evecs, evals_kinetic, P_eq, P0);
		}
		//-----Result-----
		alpha.resize(DMReals.size());
		alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * W_KM_sum * -1.0 / 9.0;
	}
}

void Lin_ESA(Matrix<double>& evecs, Vektor<double>& evals, Matrix<double>& evecs2N, Vektor<double>& evals2N, Matrix<double>& Exc2_Map, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& alpha, Vektor<double>& CIm, Vektor<double>& GReals, Vektor<double>& GImags, fftw_complex* in, fftw_complex* out, Vektor<double>& P0, Vektor<Vektor<double> >& trans_dipole_moment, Matrix< Vektor<double> >& trans_dipole_moment_2N, Vektor<double>& PM_delay){

    //-----Stimulated Emission-----

//	double PM_t = 0.0;
//	double k_21 = 0.0;
//	int low_energy = (int)evals.size() - 1;
//	int high_energy = 0;
//	int N2 = 0;
	//-----Dimer Low Energy Part-----
//	//-----Lineshape Function-----
//	DM_ESA_WZ(N2, low_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat, in, out);

//	//-----Relaxation Rate-----
//	k_21 = 2.0 * y_MNKL(low_energy, high_energy, high_energy, low_energy, evecs, parameter_vec, dipolemat)  * C_Re_WZ(evals(high_energy) - evals(low_energy) , parameter_vec(0) );
//	PM_t = P0(low_energy) + P0(high_energy) * ( 1.0 - exp(-k_21 * delay) );
//
//	//-----Result-----
//	alpha.resize(DMReals.size());
//	alpha = alpha + DMReals * inner_prod(trans_dipole_moment_2N(N2, low_energy), trans_dipole_moment_2N(N2, low_energy)) * PM_t;
//
//	//-----Dimer High Energy Part-----
//
//	//-----Lineshape Function-----
//	DM_ESA_WZ(N2, high_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat, in, out);
//
//	//-----Relaxation Rate-----
//	k_21 = 2.0 * y_MNKL(high_energy, low_energy, low_energy, high_energy, evecs, parameter_vec, dipolemat) * C_Re_WZ(evals(high_energy) - evals(low_energy) , parameter_vec(0) );
//	PM_t = P0(high_energy) * exp(-k_21 * delay) ;

	Vektor<double> DMReals;

	for(int N2 = 0 ; N2 < (int)evals2N.size(); N2++){
		for (int M = 0; M < (int)evals.size(); M++) {
			//-----Lineshape Function-----
//			DM_ESA_WZ(N2, M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat, in, out);
			DM_ESA_WZ_Markov(N2, M, DMReals,CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat);

			//-----Result-----
			alpha.resize(DMReals.size());
			alpha = alpha + DMReals * inner_prod(trans_dipole_moment_2N(N2, M), trans_dipole_moment_2N(N2, M)) * PM_delay(M) * 1.0 / 9.0;
		}
	}


//	for (int M = 0; M < (int)evals.size(); M++) {
//
//		//-----Lineshape Function-----
//		DM_WZ(M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out, p); //original 8192.0/16.0 ,8192.0 * 512.0  //Fast Version: 8192.0/256.0 ,8192.0 * 32.0 //NormalVersion 23.06.2015 8192.0/512.0 ,8192.0 * 16.0
//
//		//-----Relaxation Rate-----
//
//		if(M == low_energy){
//
//			k_21 = 2.0 * y_MNKL(M, high_energy, high_energy, M, evecs, parameter_vec, dipolemat) ;
//			PM_t = P0(low_energy) + P0(high_energy) * ( 1.0 - exp(-k_21 * delay) );
//		}
//
//		else if(M == high_energy){
//
//			k_21 = 2.0 * y_MNKL(M, low_energy, low_energy, M, evecs, parameter_vec, dipolemat);
//			PM_t = P0(high_energy) * exp(-k_21 * delay) ;
//		}
//
//
//
//		//-----Result-----
//		alpha.resize(DMReals.size());
//		alpha = alpha + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * PM_t;
//		//cout<<"M = "<<M<<"  mu-M = "<<inner_prod(trans_dipole_moment(M),trans_dipole_moment(M))<<endl;
//	}
}

double fmn_muK0(int m, int n, int K, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<Vektor<double> >& trans_dipole_moment){

	double rot_beitrag = 0;
	Vektor<double> Rmn = dipolemat(1)(m) - dipolemat(1)(n);
	Vektor<double> mu_cross = Cross(dipolemat(0)(m), dipolemat(0)(n));
	Vektor<double> mu_K = trans_dipole_moment(K);

//	rot_beitrag = 2.0 / 15.0 * inner_prod(Rmn, mu_cross) * inner_prod(mu_K, mu_K) - 1.0 / 15.0 * inner_prod(Rmn, mu_K) * inner_prod(mu_cross, mu_K);

	rot_beitrag = 1.0 / 9.0 * inner_prod(Rmn, mu_cross) * inner_prod(mu_K, mu_K);

	return rot_beitrag;
}

void CD_GSBleach(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& CD, Vektor<double>& CIm, Vektor<Vektor<double> >& trans_dipole_moment, Vektor<double>& P0, Vektor<double>& P_eq, Vektor<double>& evals_kinetic, Vektor< Vektor<double> >& ci_evecs, double delay){

    //-----CD Groundstate Bleaching-----
	Vektor<double> DMReals;

	for (int M = 0; M < (int)evals.size(); M++) {

		//-----Lineshape Function-----
		DM_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);

		double rot_strength = 0;

		for(int N = 0; N < (int)evals.size(); N++){
			for(int K = 0; K < (int)evals.size(); K++){

				for(int m = 0; m < (int)evals.size(); m++){
					for(int n = 0; n < (int)evals.size(); n++){

						rot_strength += W_KM(K, N, delay, ci_evecs, evals_kinetic, P_eq, P0, trans_dipole_moment) * evecs(M, m) * evecs(M, n) *	fmn_muK0(m, n, K, dipolemat, trans_dipole_moment);
					}
				}
			}
		}
		//-----Result-----
		CD.resize(DMReals.size());
		CD = CD + DMReals * rot_strength * -1;
	}

}

void CD_SE(Matrix<double>& evecs, Vektor<double>& evals, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& CD, Vektor<double>& CIm, Vektor<Vektor<double> >& trans_dipole_moment, Vektor<double>& P0, Vektor<double>& P_eq, Vektor<double>& evals_kinetic, Vektor< Vektor<double> >& ci_evecs, double delay){

    //-----TRCD Stimulated Emission-----
	Vektor<double> DMReals;

//	for (int M = 0; M < (int)evals.size(); M++) {
//		for(int K = 0; K < (int)evals.size(); K++){
//
//			//-----Lineshape Function-----
//			DM_Emission_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);
//
//			double rot_strength = 0;
//
//			for(int m = 0; m < (int)evals.size(); m++){
//				for(int n = 0; n < (int)evals.size(); n++){
//
//					rot_strength += evecs(M, m) * evecs(M, n) * fmn_muK0(m, n, K, dipolemat, trans_dipole_moment);
//				}
//			}
//
//			//-----Result-----
//			alpha.resize(DMReals.size());
//			alpha = alpha + DMReals * U(K) * rot_strength * -1;
//		}
//	}

	for (int M = 0; M < (int)evals.size(); M++) {

		//-----Lineshape Function-----
		DM_Emission_WZ_Markov(M, DMReals, CIm, evecs, evals, parameter_vec, dipolemat);
		double rot_strength = 0;

		for(int N = 0; N < (int)evals.size(); N++){
			for(int m = 0; m < (int)evals.size(); m++){
				for(int n = 0; n < (int)evals.size(); n++){

					rot_strength += W_KM(N, M, delay, ci_evecs, evals_kinetic, P_eq, P0, trans_dipole_moment) * evecs(M, m) * evecs(M, n) * fmn_muK0(m, n, N, dipolemat, trans_dipole_moment);
				}
			}
		}
		//-----Result-----
		CD.resize(DMReals.size());
		CD = CD + DMReals * rot_strength * -1;
	}

}

void CD_ESA(Matrix<double>& evecs, Vektor<double>& evals, Matrix<double>& evecs2N, Vektor<double>& evals2N, Matrix<double>& Exc2_Map, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor<double>& CD, Vektor<double>& CIm, Vektor<Vektor<double> >& trans_dipole_moment, Vektor<double>& P0, Vektor<double>& P_eq, Vektor<double>& evals_kinetic, Vektor< Vektor<double> >& ci_evecs, double delay){

    //-----TRCD Excited State Absorption-----
	Vektor<double> DMReals;

	for(int N2 = 0 ; N2 < (int)evals2N.size(); N2++){
		for (int M = 0; M < (int)evals.size(); M++){
			//-----Lineshape Function-----
			DM_ESA_WZ_Markov(N2, M, DMReals,CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat);

			double rot_strength = 0;

			for(int K = 0; K < (int)evals.size(); K++){

				for(int k = 0; k < (int)evals.size(); k++){
					for(int l = 0; l < (int)evals.size(); l++){

						for(int m = 0; m < (int)evals.size(); m++){
							for(int n = 0; n < (int)evals.size(); n++){

								if(k > l && m > n){

									rot_strength += W_KM(K, M, delay, ci_evecs, evals_kinetic, P_eq, P0, trans_dipole_moment) * evecs2N(N2, Exc2_Map(k,l)) * evecs2N(N2, Exc2_Map(m,n)) * evecs(M, k) * evecs(M, m) * fmn_muK0(l, n, K, dipolemat, trans_dipole_moment);
									rot_strength += W_KM(K, M, delay, ci_evecs, evals_kinetic, P_eq, P0, trans_dipole_moment) * evecs2N(N2, Exc2_Map(k,l)) * evecs2N(N2, Exc2_Map(m,n)) * 2 * evecs(M, k) * evecs(M, n) * fmn_muK0(l, m, K, dipolemat, trans_dipole_moment);
									rot_strength += W_KM(K, M, delay, ci_evecs, evals_kinetic, P_eq, P0, trans_dipole_moment) * evecs2N(N2, Exc2_Map(k,l)) * evecs2N(N2, Exc2_Map(m,n)) * evecs(M, l) * evecs(M, n) * fmn_muK0(k, m, K, dipolemat, trans_dipole_moment);
								}
							}
						}
					}
				}
			}
			//-----Result-----
			CD.resize(DMReals.size());
			CD = CD + DMReals * rot_strength * 1.0;
		}
	}

}


void CD_GSBleach_Dimerold(Matrix<double> evecs, Vektor<double> evals, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor<double>& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p, Vektor<double> P0, Vektor<Vektor<double> > trans_dipole_moment, Vektor<Vektor<double> > magn_dipole_moment){
	//-----Pump-Probe-CD-----
    //-----Groundstate Bleaching-----

	Vektor<double> DMReals;
	Vektor<double> CrossResult(3);

	int low_energy = (int)evals.size() - 1;
	int high_energy = 0;

	double GB_LE_1 = 0;
	double GB_HE_1 = 0;
	double GB_cross_LE = 0;
	double GB_cross_HE = 0;

	//-----Dimer Low Energy Part-----

//	//	//-----DEBUG-----
//	cout <<"DEBUG -- GSB LE Part " << endl;

	//-----4/15 (m1 . mu1) (mu1. mu1)-----
	GB_LE_1 = 4.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy) , trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment(low_energy), trans_dipole_moment(low_energy)) * P0(low_energy);

//	cout <<"DEBUG -- LE_1: " << GB_LE_1 << " P0: "<< P0(low_energy) << endl;

	//-----3/15 (m1 . mu1) (mu2. mu2) + 1/15 (m1 . mu2) (mu1 . mu2)-----
	GB_cross_LE = (3.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy) , trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) +
				   1.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy) , trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment(low_energy), trans_dipole_moment(high_energy)) ) * P0(high_energy);

//	cout <<"DEBUG -- Cross_LE: " << GB_cross_LE << " P0: " << P0(high_energy) << endl;

	DM_WZ(low_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);

	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * (GB_LE_1 + GB_cross_LE);

	//-----Dimer High Energy Part-----

	//-----4/15 (m2 . mu2) (mu2. mu2)-----
	GB_HE_1 = 4.0 / 15.0 * inner_prod(magn_dipole_moment(high_energy) , trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) * P0(high_energy);

//	cout <<"DEBUG -- HE_1: " << GB_HE_1 << endl;


	//-----3/15 (m2 . mu2) (mu1. mu1) + 1/15 (m2 . mu1) (mu2 . mu1)-----
	GB_cross_HE = (3.0 / 15.0 * inner_prod(magn_dipole_moment(high_energy) , trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment(low_energy), trans_dipole_moment(low_energy)) +
			   	   1.0 / 15.0 * inner_prod(magn_dipole_moment(high_energy) , trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(low_energy)) ) * P0(low_energy);

//	cout <<"DEBUG -- Cross_HE: " << GB_cross_HE << endl;

	DM_WZ(high_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);

	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * (GB_HE_1 + GB_cross_HE);

}

void CD_StimEmission_Dimerold(Matrix<double> evecs, Vektor<double> evals, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor<double>& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p, Vektor<double> P0, Vektor<Vektor<double> > trans_dipole_moment, Vektor<Vektor<double> > magn_dipole_moment, double delay){
	//-----Pump-Probe-CD-----
    //-----Groundstate Bleaching-----

	Vektor<double> DMReals;

	int low_energy = (int)evals.size() - 1;
	int high_energy = 0;

	double k_21 = 0;
	double PM_t = 0;

	double GB_LE_1 = 0;
	double GB_HE_1 = 0;
	double GB_cross_LE = 0;
	double GB_cross_HE = 0;

	//-----Dimer Low Energy Part-----

//	//	//-----DEBUG-----
//	cout <<"DEBUG -- SE LE Part " << endl;

	//-----Lineshape Function-----
	DM_Emission_WZ(low_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);

	//-----Relaxation Rate-----
	k_21 = 2.0 * y_MNKL(low_energy, high_energy, high_energy, low_energy, evecs, parameter_vec, dipolemat) * C_Re_WZ(evals(high_energy) - evals(low_energy), parameter_vec(0));
	PM_t = (1 - exp(-k_21 * delay)) * P0(high_energy);

	//-----4/15 (m1 . mu1) (mu1. mu1) U(w1)-----
	GB_LE_1 = 4.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy) , trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment(low_energy), trans_dipole_moment(low_energy)) * P0(low_energy);

	//-----( 3/15 (m1 . mu1) (mu2. mu2) + 1/15 (m1 . mu2) (mu1 . mu2) ) U(w2) (1-Exp[-k21 t])-----
	GB_cross_LE = (3.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy) , trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) +
				   1.0 / 15.0 * inner_prod(magn_dipole_moment(low_energy), trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(low_energy)) ) * PM_t;

//	//	//-----DEBUG-----
//	cout << "DEBUG -- LE Part: " << GB_LE_1 << endl;
//	cout << "DEBUG -- LE Part: " << GB_cross_LE <<" PM(t): " << PM_t << endl;

	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * (GB_LE_1 + GB_cross_LE);

	//-----Dimer High Energy Part-----

//	//	//-----DEBUG-----
//	cout <<"DEBUG -- SE HE Part " << endl;

	//-----Lineshape Function-----
	DM_Emission_WZ(high_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);

	//-----Relaxation Rate-----
	k_21 = 2.0 * y_MNKL(high_energy, low_energy, low_energy, high_energy, evecs, parameter_vec, dipolemat) * C_Re_WZ(evals(high_energy) - evals(low_energy), parameter_vec(0));
	PM_t = P0(high_energy) * exp(-k_21 * delay);

	//-----4/15 (m2 . mu2) (mu2. mu2) U(w2) Exp[-k21 t]-----
	GB_HE_1 = 4.0 / 15.0 * inner_prod(magn_dipole_moment(high_energy) , trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) * PM_t;

//	//	//-----DEBUG-----
//	cout << "DEBUG -- HE Part: " << GB_HE_1 <<" PM(t): " << PM_t << endl;

	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * (GB_HE_1);

}

void CD_ESA_Dimerold(Matrix<double> evecs, Vektor<double> evals, Matrix<double> evecs2N, Vektor<double> evals2N, Matrix<double> Exc2_Map, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor<double>& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p, Vektor<double> P0, Vektor<Vektor<double> > trans_dipole_moment, Matrix< Vektor<double> > trans_dipole_moment_2N, Vektor<Vektor<double> > magn_dipole_moment, Matrix< Vektor<double> > magn_dipole_moment_2N, double delay){

    //-----Excited State Absorption-----

	Vektor<double> DMReals;
	double PM_t = 0.0;
	double k_21 = 0.0;
	int low_energy = (int)evals.size() - 1;
	int high_energy = 0;
	int N2 = 0;

	double GB_cross_HE = 0;
	double GB_cross_LE_1 = 0;
	double GB_cross_LE_2 = 0;

//	//	-----DEBUG-----
//	cout <<"DEBUG -- ESA LE Part " << endl;

	//-----Dimer Low Energy Part-----

	//-----Lineshape Function-----
	DM_ESA_WZ(N2, low_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat, in, out);

	//-----Relaxation Rate-----
	k_21 = 2.0 * y_MNKL(low_energy, high_energy, high_energy, low_energy, evecs, parameter_vec, dipolemat)  * C_Re_WZ(evals(high_energy) - evals(low_energy) , parameter_vec(0) );
	PM_t = P0(high_energy) * ( 1.0 - exp(-k_21 * delay) );

	//-----( 3/15 (m1,2N . mu1,2N) (mu2. mu2) + 1/15 (m1,2N . mu2) (mu1,2N . mu2) ) U(w2) (1-Exp[-k21 t])-----
	GB_cross_LE_1 = (3.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, low_energy), trans_dipole_moment_2N(N2, low_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) +
			         1.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, low_energy), trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment_2N(N2, low_energy), trans_dipole_moment(high_energy))) * PM_t;


	//-----( 3/15 (m1,2N . mu1,2N) (mu2. mu2) + 1/15 (m1,2N . mu2) (mu1,2N . mu2) ) U(w2) (1-Exp[-k21 t])-----
	GB_cross_LE_2 = (3.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, low_energy), trans_dipole_moment_2N(N2, low_energy)) * inner_prod(trans_dipole_moment(low_energy), trans_dipole_moment(low_energy)) +
			         1.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, low_energy), trans_dipole_moment(low_energy)) * inner_prod(trans_dipole_moment_2N(N2, low_energy), trans_dipole_moment(low_energy))) * P0(low_energy);


//	//	//-----DEBUG-----
//	cout << "DEBUG -- LE Part: " << GB_cross_LE_1 << " PM(t): " << PM_t << endl;
//	cout << "DEBUG -- LE Part: " << GB_cross_LE_2 << " PM(0): " << P0 << endl;

	//-----Result-----
	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * (GB_cross_LE_1 + GB_cross_LE_2);


	//-----Dimer High Energy Part-----

	//-----Lineshape Function-----
	DM_ESA_WZ(N2, high_energy, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat, in, out);

	//-----Relaxation Rate-----
	k_21 = 2.0 * y_MNKL(high_energy, low_energy, low_energy, high_energy, evecs, parameter_vec, dipolemat) * C_Re_WZ(evals(high_energy) - evals(low_energy) , parameter_vec(0) );
	PM_t = P0(high_energy) * exp(-k_21 * delay) ;

	//-----( 3/15 (m2,2N . mu2,2N) (mu2. mu2) + 1/15 (m2,2N . mu2) (mu2,2N . mu2) ) U(w2) (1-Exp[-k21 t])-----
	GB_cross_HE = (3.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, high_energy), trans_dipole_moment_2N(N2, high_energy)) * inner_prod(trans_dipole_moment(high_energy), trans_dipole_moment(high_energy)) +
				   1.0 / 15.0 * inner_prod(magn_dipole_moment_2N(N2, high_energy), trans_dipole_moment(high_energy)) * inner_prod(trans_dipole_moment_2N(N2, high_energy), trans_dipole_moment(high_energy))) * PM_t;

//	cout << "DEBUG -- HE Part: " << GB_cross_HE << " PM(t): " << PM_t << endl;

	//-----Result-----
	alpha.resize(DMReals.size());
	alpha = alpha + DMReals * GB_cross_HE;

}

void ACD(Matrix<double>& evecs, Vektor<double>& evals, Vektor< ChlAtoms >& Chl_Vec, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& parameter_vec, Vektor< Vektor<double> >& ACD_Edge, Vektor< Vektor<double> >& ACD_Face, Vektor< Vektor<double> >& CD, Vektor< Vektor<double> >& alpha, Vektor<double>& CIm, Vektor<double>& GReals, Vektor<double>& GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p){

    //Allocate a RM Vectors for all RM values of eigenvalues M
	Vektor<double> RM(evals.size());
	for(int i = 0; i < (int) RM.size(); i++){
		RM(i) = 0;
	}

	Vektor<double> RM_Edge(evals.size());
	for(int i = 0; i < (int) RM_Edge.size(); i++){
		RM_Edge(i) = 0;
	}

	Vektor<double> RM_Face(evals.size());
	for(int i = 0; i < (int) RM_Face.size(); i++){
		RM_Face(i) = 0;
	}

    Vektor<double> CrossResult(3);
//	int LE = (int)evals.size() - 1;
//	int HE = 0;

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}

//	cout << trans_dipole_moment << endl;
    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size(); M++){

		for (int i = 0; i < (int) evals.size() - 1; i++) {
			for (int j = i + 1; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);

				RM(M) += 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

				//-----R-Dependent Version!-----
//				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

//				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);
//				RM_Edge(M) += 1.0/2.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
			}
		}

		for (int i = 0; i < (int) evals.size(); i++) {
			for (int j = 0; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);

//				RM(M) += 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);

//				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

//				RM_Edge(M) += 1.0/2.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
			}
		}

//		cout << "RM(" << M <<") = " << RM(M) << endl;
//		cout << "Centers: " << Chl_Vec(0).Mol_Center << "  " << Chl_Vec(1).Mol_Center << endl;
//		cout << "Dip Vec: " << Chl_Vec(0).Dipole_Mom << "  " << Chl_Vec(1).Dipole_Mom << endl;
//
//		cout << "M = " << M << " RM = " << RM(M) << endl;
//
//		cout << "Quadrupole Moment 1:" << Chl_Vec(0).Quadrupole_Moment << endl;
//		cout << "Quadrupole Moment 2:" << Chl_Vec(1).Quadrupole_Moment << endl;

		//-----The ACD Part-----
    	for(int i = 0; i < (int) evals.size(); i++){
    		for(int j = 0; j < (int) evals.size(); j++){

    			//-----R-Dependent Version!-----
//    			crossprod(Chl_Vec(j).Mol_Center, Chl_Vec(j).Dipole_Mom, CrossResult);
//       		RM_Face(M) +=  -1.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
//       		RM_Face(M) +=   1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

    			RM_Face(M) +=   1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(1,2) );

//    			RM_Edge(M) +=   1.0/2.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
//    			RM_Edge(M) +=  -1.0/2.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

//    			RM_Face(M) +=  1.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
//    			RM_Edge(M) += -1.0 / 2.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

//    			cout << "ACD PART M :" << M << endl;
//
//    			cout << evecs(M, i) << "  "  << evecs(M, j) << "  "  << Chl_Vec(i).Dipole_Mom(1) << "  "  << Chl_Vec(j).Quadrupole_Moment(0,2) << endl;
//    			cout << evecs(M, i)  << "  "  <<  evecs(M, j)  << "  "  <<  Chl_Vec(i).Dipole_Mom(0)  << "  "  <<  Chl_Vec(j).Quadrupole_Moment(1,2) << endl;

    		}
    	}

//    	//-----Isotropic part of rM_edge-----
//
//		for (int i = 0; i < (int) evals.size() - 1; i++) {
//			for (int j = i + 1; j < (int) evals.size(); j++) {
//
//				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);
//
//				RM(M) += 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
//
//				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
//				RM_Edge(M) += 1.0/2.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
//			}
//		}
//
////		cout << "RM(" << M <<") = " << RM(M) << endl;
////		cout << "Centers: " << Chl_Vec(0).Mol_Center << "  " << Chl_Vec(1).Mol_Center << endl;
////		cout << "Dip Vec: " << Chl_Vec(0).Dipole_Mom << "  " << Chl_Vec(1).Dipole_Mom << endl;
////
////		cout << "M = " << M << " RM = " << RM(M) << endl;
////
////		cout << "Quadrupole Moment 1:" << Chl_Vec(0).Quadrupole_Moment << endl;
////		cout << "Quadrupole Moment 2:" << Chl_Vec(1).Quadrupole_Moment << endl;
//
//		//-----The ACD Part-----
//    	for(int i = 0; i < (int) evals.size(); i++){
//    		for(int j = 0; j < (int) evals.size(); j++){
//
//    			crossprod(Chl_Vec(j).Mol_Center, Chl_Vec(j).Dipole_Mom, CrossResult);
//
//    			RM_Face(M) +=  -1.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
//    	    	RM_Face(M) +=   1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
//
//    			RM_Edge(M) +=   1.0/2.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
//    			RM_Edge(M) +=  -1.0/2.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
//
////    			RM_Face(M) +=  1.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
////    			RM_Edge(M) += -1.0 / 2.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
//
////    			cout << "ACD PART M :" << M << endl;
////
////    			cout << evecs(M, i) << "  "  << evecs(M, j) << "  "  << Chl_Vec(i).Dipole_Mom(1) << "  "  << Chl_Vec(j).Quadrupole_Moment(0,2) << endl;
////    			cout << evecs(M, i)  << "  "  <<  evecs(M, j)  << "  "  <<  Chl_Vec(i).Dipole_Mom(0)  << "  "  <<  Chl_Vec(j).Quadrupole_Moment(1,2) << endl;
//
//    		}
//    	}

//    	cout << "RM_Face(" << M << ") = " << RM_Face(M) << endl;
//    	cout << "RM_Edge(" << M << ") = " << RM_Edge(M) << endl;
//    	cout << "Evecs: " << evecs(M, 0) << "  " << evecs(M, 1) << endl;

        //CHANGE LINESHAPEFUNCTION DM_WZ HERE!
		DM_WZ(M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out);

        //-----CD Part-----
        CD(0).resize(DMReals.size());
        CD(0) = CD(0) + DMReals * RM(M);

        //-----ACD Edge & Face Part-----
        ACD_Edge(0).resize(DMReals.size());
        ACD_Edge(0) = ACD_Edge(0) + DMReals * RM_Edge(M);

        ACD_Face(0).resize(DMReals.size());
        ACD_Face(0) = ACD_Face(0) + DMReals * RM_Face(M);

        //-----ABSORPTION PART-----
        alpha(0).resize(DMReals.size());
        alpha(0) = alpha(0) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));

//        if(M == LE){
//			//-----Low Energy Contribution-----
//			CD(1).resize(DMReals.size());
//			CD(1) = CD(1) + DMReals * RM(M);
//
//			ACD_Edge(1).resize(DMReals.size());
//			ACD_Edge(1) = ACD_Edge(1) + DMReals * RM_Edge(M);
//
//			ACD_Face(1).resize(DMReals.size());
//			ACD_Face(1) = ACD_Face(1) + DMReals * RM_Face(M);
//
//			alpha(1).resize(DMReals.size());
//			alpha(1) = alpha(1) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
//        }
//
//        else if(M == HE){
//			//-----High Energy Contribution-----
//			CD(2).resize(DMReals.size());
//			CD(2) = CD(2) + DMReals * RM(M);
//
//			ACD_Edge(2).resize(DMReals.size());
//			ACD_Edge(2) = ACD_Edge(2) + DMReals * RM_Edge(M);
//
//			ACD_Face(2).resize(DMReals.size());
//			ACD_Face(2) = ACD_Face(2) + DMReals * RM_Face(M);
//
//			alpha(2).resize(DMReals.size());
//			alpha(2) = alpha(2) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
//        }

    }
}

void get_ACD_RotStrength(Matrix<double> evecs, Vektor<double> evals, Vektor< ChlAtoms > Chl_Vec, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor< Vektor<double> >& RotStrength_vec){

    //Allocate a RM Vectors for all RM values of eigenvalues M
	Vektor<double> RM(evals.size());
	for(int i = 0; i < (int) RM.size(); i++){
		RM(i) = 0;
	}

	//-----RM Edge & Face-----
	Vektor<double> RM_Edge(evals.size());
	for(int i = 0; i < (int) RM_Edge.size(); i++){
		RM_Edge(i) = 0;
	}

	Vektor<double> RM_Face(evals.size());
	for(int i = 0; i < (int) RM_Face.size(); i++){
		RM_Face(i) = 0;
	}
	//-----RM Edge & Face z-Part-----
	Vektor<double> RM_zPart(evals.size());
	for(int i = 0; i < (int) RM_zPart.size(); i++){
		RM_zPart(i) = 0;
	}

	//-----RM Edge & Face Quadrupole Part-----
	Vektor<double> RM_Quadrupole(evals.size());
	for(int i = 0; i < (int) RM_Quadrupole.size(); i++){
		RM_Quadrupole(i) = 0;
	}

	//-----RM Edge & Face Quadrupole Part-----
	Vektor<double> RM_Diff_Part1(evals.size());
	for(int i = 0; i < (int) RM_Diff_Part1.size(); i++){
		RM_Diff_Part1(i) = 0;
	}

	//-----RM Edge & Face Quadrupole Part-----
	Vektor<double> RM_Diff_Part2(evals.size());
	for(int i = 0; i < (int) RM_Diff_Part2.size(); i++){
		RM_Diff_Part2(i) = 0;
	}

	//-----RM Edge & Face Quadrupole Part-----
	Vektor<double> RM_Diff_Part3(evals.size());
	for(int i = 0; i < (int) RM_Diff_Part3.size(); i++){
		RM_Diff_Part3(i) = 0;
	}

	//-----RM Edge & Face Intrinsic Quadrupole Part-----
	Vektor<double> RM_Intrinsic_Quadrupole(evals.size());
	for(int i = 0; i < (int) RM_Intrinsic_Quadrupole.size(); i++){
		RM_Intrinsic_Quadrupole(i) = 0;
	}

    Vektor<double> CrossResult(3);
	int LE = (int)evals.size() - 1;
	int HE = 0;

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}

//    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size(); M++){

    	//-----Isotropic part of rM_edge-----

		for (int i = 0; i < (int) evals.size() - 1; i++) {
			for (int j = i + 1; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);

				RM(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
				RM_Edge(M) += 1.0/2.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
			}
		}

//		cout << "RM(" << M <<") = " << RM(M) << endl;
//		cout << "Centers: " << Chl_Vec(0).Mol_Center << "  " << Chl_Vec(1).Mol_Center << endl;
//		cout << "Dip Vec: " << Chl_Vec(0).Dipole_Mom << "  " << Chl_Vec(1).Dipole_Mom << endl;
//
//		cout << "M = " << M << " RM = " << RM(M) << endl;
//
//		cout << "Quadrupole Moment 1:" << Chl_Vec(0).Quadrupole_Moment << endl;
//		cout << "Quadrupole Moment 2:" << Chl_Vec(1).Quadrupole_Moment << endl;

		//-----The ACD Part-----
    	for(int i = 0; i < (int) evals.size(); i++){
    		for(int j = 0; j < (int) evals.size(); j++){

    			crossprod(Chl_Vec(j).Mol_Center, Chl_Vec(j).Dipole_Mom, CrossResult);

    			RM_Face(M) +=  -1.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
    	    	RM_Face(M) +=   1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

    			RM_Edge(M) +=   1.0/2.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
    			RM_Edge(M) +=  -1.0/2.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

    			//-----ACD Rotational-Strength Parts-----
    			RM_zPart(M) +=  1.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * CrossResult(2);
    			RM_Quadrupole(M) += 1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
    			RM_Intrinsic_Quadrupole(M) += 1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(1,2) );

//    			crossprod( (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center) , Chl_Vec(j).Dipole_Mom, CrossResult);
    			RM_Diff_Part1(M) +=  1.0 / 2.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(2) * Cross((Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center), Chl_Vec(j).Dipole_Mom)(2);
    			RM_Diff_Part2(M) += -1.0 / 2.0 * evecs(M, i) * evecs(M, j) * Chl_Vec(j).Dipole_Mom(2) * Cross((Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center), Chl_Vec(i).Dipole_Mom)(2);
    			RM_Diff_Part3(M) +=  1.0 / 2.0 * evecs(M, i) * evecs(M, j) * (Chl_Vec(j).Mol_Center - Chl_Vec(i).Mol_Center)(2) * Cross(Chl_Vec(j).Dipole_Mom, Chl_Vec(i).Dipole_Mom)(2);


//    			RM_Face(M) +=  1.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
//    			RM_Edge(M) += -1.0 / 2.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );

//    			cout << "ACD PART M :" << M << endl;
//
//    			cout << evecs(M, i) << "  "  << evecs(M, j) << "  "  << Chl_Vec(i).Dipole_Mom(1) << "  "  << Chl_Vec(j).Quadrupole_Moment(0,2) << endl;
//    			cout << evecs(M, i)  << "  "  <<  evecs(M, j)  << "  "  <<  Chl_Vec(i).Dipole_Mom(0)  << "  "  <<  Chl_Vec(j).Quadrupole_Moment(1,2) << endl;

    		}
    	}

//    	cout << "RM(M) : " << RM(M) << endl << "RM z-Part(M) : " << RM_zPart(M) << endl << "RM Quadrupole(M) : " << RM_Quadrupole(M) << endl << "RM_Intrinsic Quadrupole(M) : " << RM_Intrinsic_Quadrupole(M) << endl;

//    	cout << "RM_Face(" << M << ") = " << RM_Face(M) << endl;
//    	cout << "RM_Edge(" << M << ") = " << RM_Edge(M) << endl;
//    	cout << "Evecs: " << evecs(M, 0) << "  " << evecs(M, 1) << endl;

    	RotStrength_vec.resize(9);
    	for(int i = 0; i < (int)RotStrength_vec.size(); i++)
    		RotStrength_vec(i).resize( (int)RM.size() );

        if(M == LE){
			//-----Low Energy Contribution-----
        	RotStrength_vec(0)(LE) = RM(M);
        	RotStrength_vec(1)(LE) = RM_zPart(M);
        	RotStrength_vec(2)(LE) = RM_Quadrupole(M);
        	RotStrength_vec(3)(LE) = RM_Intrinsic_Quadrupole(M);
        	RotStrength_vec(4)(LE) = RM(M) - RM_zPart(M) + RM_Quadrupole(M);
        	RotStrength_vec(5)(LE) = RM(M) + RM_zPart(M) - RM_Quadrupole(M);

        	RotStrength_vec(6)(LE) = RM_Diff_Part1(M);
        	RotStrength_vec(7)(LE) = RM_Diff_Part2(M);
        	RotStrength_vec(8)(LE) = RM_Diff_Part3(M);

        }

        else if(M == HE){
			//-----High Energy Contribution-----
        	RotStrength_vec(0)(HE) = RM(M);
        	RotStrength_vec(1)(HE) = RM_zPart(M);
        	RotStrength_vec(2)(HE) = RM_Quadrupole(M);
        	RotStrength_vec(3)(HE) = RM_Intrinsic_Quadrupole(M);
        	RotStrength_vec(4)(HE) = RM(M) - RM_zPart(M) + RM_Quadrupole(M);
        	RotStrength_vec(5)(HE) = RM(M) + RM_zPart(M) - RM_Quadrupole(M);

          	RotStrength_vec(6)(HE) = RM_Diff_Part1(M);
          	RotStrength_vec(7)(HE) = RM_Diff_Part2(M);
          	RotStrength_vec(8)(HE) = RM_Diff_Part3(M);
        }

    }

}

void ACD_face(Matrix<double> evecs, Vektor<double> evals, Vektor< ChlAtoms > Chl_Vec, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor< Vektor<double> >& CD, Vektor< Vektor<double> >& alpha, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p){

    //Allocate a RM Vector for all RM values of eigenvalues M
	Vektor<double> RM(evals.size());
	for(int i = 0; i < (int) RM.size(); i++){
		RM(i) = 0;
	}

    Vektor<double> CrossResult(3);
	int LE = (int)evals.size() - 1;
	int HE = 0;

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}

    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size() ;M++){

    	//-----Isotropic part of rM_edge-----

		for (int i = 0; i < (int) evals.size() - 1; i++) {
			for (int j = i + 1; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);

				RM(M) += 1.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
			}
		}

		//-----The ACD_edge Part-----
    	for(int i = 0; i < (int) evals.size() - 1; i++){
    		for(int j = 0; j < (int) evals.size() - 1; j++){

    			//----- (R_k x mu_k)(z) * mu_l(z)-----
    			RM(M) += 1.0 * ( evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Quadrupole_Moment(0,2) - evecs(M, i) * evecs(M, j) * Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Quadrupole_Moment(1,2) );
    		}
    	}


        //CHANGE LINESHAPEFUNCTION DM_WZ HERE!
        DM_WZ(M,DMReals,GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,CIm,evecs,evals,parameter_vec,dipolemat, in, out);

        //-----CD Part-----
        CD(0).resize(DMReals.size());
        CD(0) = CD(0) + DMReals * RM(M);

        //-----ABSORPTION PART-----
        alpha(0).resize(DMReals.size());
        alpha(0) = alpha(0) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));

        if(M == LE){
			//-----Low Energy Contribution-----
			CD(1).resize(DMReals.size());
			CD(1) = CD(1) + DMReals * RM(M);

			alpha(1).resize(DMReals.size());
			alpha(1) = alpha(1) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
        }

        else if(M == HE){
			//-----High Energy Contribution-----
			CD(2).resize(DMReals.size());
			CD(2) = CD(2) + DMReals * RM(M);

			alpha(2).resize(DMReals.size());
			alpha(2) = alpha(2) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M));
        }

    }
}

void Linear_Dichroism(Matrix<double> evecs, Vektor<double> evals, Vektor< ChlAtoms > Chl_Vec, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor<double> parameter_vec, Vektor< Vektor<double> >& LD, Vektor<double> CIm, Vektor<double> GReals, Vektor<double> GImags, fftw_complex* in, fftw_complex* out, fftw_plan& p){

	Vektor<double> Symmetrie_Achse(3);

	Symmetrie_Achse(0) = 0.0;
	Symmetrie_Achse(1) = 0.0;
	Symmetrie_Achse(2) = 1.0;

    Vektor<double> CrossResult(3);
//	int LE = (int)evals.size() - 1;
//	int HE = 0;

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}

    Vektor<double> DMReals;
    double Winkel = 0.0;

    for(int M = 0; M < (int)trans_dipole_moment.size(); M++){

    	//-----Winkel zwischen Exc. Dipolmoment & Symmetrieachse-----

    	Winkel = acos(inner_prod( trans_dipole_moment(M), Symmetrie_Achse ) / ( norm_2(trans_dipole_moment(M)) ));

    	if(norm_2(trans_dipole_moment(M)) < 1e-12){
    		continue;
    	}

//    	cout << "RM_Edge(" << M << ") = " << RM_Edge(M) << endl;
//    	cout << "Evecs: " << evecs(M, 0) << "  " << evecs(M, 1) << endl;

//        //CHANGE LINESHAPEFUNCTION DM_WZ HERE!
		DM_WZ(M, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs, evals, parameter_vec, dipolemat, in, out); //original 8192.0/16.0 ,8192.0 * 512.0  //Fast Version: 8192.0/256.0 ,8192.0 * 32.0 //NormalVersion 23.06.2015 8192.0/512.0 ,8192.0 * 16.0

		//-----LD Part-----
		LD(0).resize(DMReals.size());
		LD(0) = LD(0) + DMReals * inner_prod(trans_dipole_moment(M),trans_dipole_moment(M)) * (1 - 3 * cos(Winkel) * cos(Winkel) );

//		if (M == LE) {
//			//-----Low Energy Contribution-----
//			LD(1).resize(DMReals.size());
//			LD(1) = LD(1) + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * (1 - 3 * cos(Winkel) * cos(Winkel));
//		}
//
//		else if (M == HE) {
//			//-----High Energy Contribution-----
//			LD(2).resize(DMReals.size());
//			LD(2) = LD(2) + DMReals * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * (1 - 3 * cos(Winkel) * cos(Winkel));
//		}

    }

}

void get_LD_RotStrength(Matrix<double> evecs, Vektor<double> evals, Vektor< ChlAtoms > Chl_Vec, Vektor< Vektor< Vektor<double> > > dipolemat, Vektor< Vektor<double> >& LDStrength_vec){

    //Allocate a RM Vectors for all RM values of eigenvalues M
	Vektor<double> RM(evals.size());
	for(int i = 0; i < (int) RM.size(); i++){
		RM(i) = 0;
	}

    Vektor<double> CrossResult(3);
	int LE = (int)evals.size() - 1;
	int HE = 0;

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}


	Vektor<double> Symmetrie_Achse(3);

	Symmetrie_Achse(0) = 0.0;
	Symmetrie_Achse(1) = 0.0;
	Symmetrie_Achse(2) = 1.0;

	double Winkel = 0.0;

//    Vektor<double> DMReals;

    for(int M = 0; M < (int)RM.size(); M++){

    	//-----Winkel zwischen Exc. Dipolmoment & Symmetrieachse-----
    	Winkel = acos(inner_prod(trans_dipole_moment(M), Symmetrie_Achse) / (norm_2(trans_dipole_moment(M))));

    	LDStrength_vec.resize(9);
    	for(int i = 0; i < (int)LDStrength_vec.size(); i++)
    		LDStrength_vec(i).resize( (int)RM.size() );

        if(M == LE){
			//-----Low Energy Contribution-----
        	LDStrength_vec(0)(LE) = inner_prod(trans_dipole_moment(M),trans_dipole_moment(M)) * (1 - 3 * cos(Winkel) * cos(Winkel) );

//        	RotStrength_vec(0)(LE) = RM(M);
//        	RotStrength_vec(1)(LE) = RM_zPart(M);
//        	RotStrength_vec(2)(LE) = RM_Quadrupole(M);
//        	RotStrength_vec(3)(LE) = RM_Intrinsic_Quadrupole(M);
//        	RotStrength_vec(4)(LE) = RM(M) - RM_zPart(M) + RM_Quadrupole(M);
//        	RotStrength_vec(5)(LE) = RM(M) + RM_zPart(M) - RM_Quadrupole(M);
//
//        	RotStrength_vec(6)(LE) = RM_Diff_Part1(M);
//        	RotStrength_vec(7)(LE) = RM_Diff_Part2(M);
//        	RotStrength_vec(8)(LE) = RM_Diff_Part3(M);

        }

        else if(M == HE){
			//-----High Energy Contribution-----
        	LDStrength_vec(0)(HE) = inner_prod(trans_dipole_moment(M),trans_dipole_moment(M)) * (1 - 3 * cos(Winkel) * cos(Winkel) );

//        	RotStrength_vec(0)(HE) = RM(M);
//        	RotStrength_vec(1)(HE) = RM_zPart(M);
//        	RotStrength_vec(2)(HE) = RM_Quadrupole(M);
//        	RotStrength_vec(3)(HE) = RM_Intrinsic_Quadrupole(M);
//        	RotStrength_vec(4)(HE) = RM(M) - RM_zPart(M) + RM_Quadrupole(M);
//        	RotStrength_vec(5)(HE) = RM(M) + RM_zPart(M) - RM_Quadrupole(M);
//
//          	RotStrength_vec(6)(HE) = RM_Diff_Part1(M);
//          	RotStrength_vec(7)(HE) = RM_Diff_Part2(M);
//          	RotStrength_vec(8)(HE) = RM_Diff_Part3(M);
        }

    }
}

void get_RotStrengths(Matrix<double>& evecs, Vektor<double>& evals, Vektor< ChlAtoms >& Chl_Vec, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& LDStrength_vec, Vektor<double>& RotStrength_vec, Vektor<double>& ACD_RotStrength_vec, Vektor<double>& CIm, Vektor<double>& parameter_vec, Vektor<double>& wM0Strength_vec, Vektor < Vektor <double> >& Rot_Beitrag, Vektor < Vektor <double> >& ACDRot_Beitrag){

    //-----Allocate Strength-Vectors for all eigenvalues M-----

	Vektor<double> RM(evals.size());
	for(int i = 0; i < (int) RM.size(); i++){
		RM(i) = 0;
	}

	Vektor<double> RM_Face(evals.size());
	for(int i = 0; i < (int) RM_Face.size(); i++){
		RM_Face(i) = 0;
	}

	Vektor<double> LD_Strength(evals.size());
	for(int i = 0; i < (int) RM_Face.size(); i++){
		LD_Strength(i) = 0;
	}

    Vektor<double> CrossResult(3);

	//ABSORPTION-PART
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);
	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
	}

	for (int i = 0; i < (int) trans_dipole_moment.size(); i++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {
			trans_dipole_moment(i) += dipolemat(0)(j) * evecs(i, j);
		}
	}

	Vektor<double> wM0_vec(evals.size());
	for(int i = 0; i < (int) wM0_vec.size(); i++){
		wM0_vec(i) = 0;
	}

//	Vektor < Vektor <double> >  Rot_Beitrags(3);



    for(int M = 0; M < (int)RM.size(); M++){

		for (int i = 0; i < (int) evals.size() - 1; i++) {
			for (int j = i + 1; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);

				RM(M) += 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

				//-----Pick a single (lowest?) Exciton-State-----
				if(M == 21){

					Vektor<double> tmp(6);
					tmp(0) = i;
					tmp(1) = j;
					tmp(2) = evecs(M, i);
					tmp(3) = evecs(M, j);
					tmp(4) = inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);
					tmp(5) = 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult);

					push_back(Rot_Beitrag, tmp);

//					output << i << "  " << j << "  " << evecs(M, i) << "  " << evecs(M, j) << "  " << Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center << "  " << CrossResult;
//					output << "  " << inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult) << "  " << 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult) << endl;
				}

			}
		}

		for (int i = 0; i < (int) evals.size(); i++) {
			for (int j = 0; j < (int) evals.size(); j++) {

				crossprod(Chl_Vec(i).Dipole_Mom, Chl_Vec(j).Dipole_Mom, CrossResult);
				RM_Face(M) += 1.0 * evecs(M, i) * evecs(M, j) * (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);

				if (M == 21) {

					Vektor<double> tmp(7);
					tmp(0) = i;
					tmp(1) = j;
					tmp(2) = evecs(M, i);
					tmp(3) = evecs(M, j);
					tmp(4) = (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);
					tmp(5) = evecs(M, i) * evecs(M, j) * (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);
					tmp(6) = evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(1,2) );

					push_back(ACDRot_Beitrag, tmp);

					//					output << i << "  " << j << "  " << evecs(M, i) << "  " << evecs(M, j) << "  " << Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center << "  " << CrossResult;
					//					output << "  " << inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult) << "  " << 2.0 / 3.0 * evecs(M, i) * evecs(M, j) * inner_prod(Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center, CrossResult) << endl;
				}

			}
		}

		//-----The ACD Part-----
    	for(int i = 0; i < (int) evals.size(); i++){
    		for(int j = 0; j < (int) evals.size(); j++){

    			RM_Face(M) +=   1.0 * evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(1,2) );

//    			if (M == 23) {
//
//					tmp(6) = evecs(M, i) * evecs(M, j) * (Chl_Vec(i).Mol_Center - Chl_Vec(j).Mol_Center)(2) * CrossResult(2);
//					tmp(7) = evecs(M, i) * evecs(M, j) * ( Chl_Vec(i).Dipole_Mom(1) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(0,2) - Chl_Vec(i).Dipole_Mom(0) * Chl_Vec(j).Intrinsic_Quadrupole_Moment(1,2) );
//					push_back(Rot_Beitrag, tmp);
//    			}

    		}
    	}

//    	double wM0 = wM0_WZ(M, evals, evecs, parameter_vec, CIm, dipolemat);

//    	wM0_vec(M) = wM0;

//    	cout <<"evals: " << M << "  " << evals(M) << "  wM0: " << wM0 << "  " << wM0_vec(M) << endl;
//    	cout << RM(M) << "  " << RM_Face(M) << endl << endl;

    }

    Vektor<double> Symmetrie_Achse(3);

	Symmetrie_Achse(0) = 0.0;
	Symmetrie_Achse(1) = 0.0;
	Symmetrie_Achse(2) = 1.0;

	double Winkel = 0.0;

	for (int M = 0; M < (int) RM.size(); M++) {

		//-----Winkel zwischen Exc. Dipolmoment & Symmetrieachse-----
		Winkel = acos(inner_prod(trans_dipole_moment(M), Symmetrie_Achse) / (norm_2(trans_dipole_moment(M))));

		LD_Strength(M) += inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * ( 1 - 3 * cos(Winkel) * cos(Winkel) );

	}

	RotStrength_vec = RM;
	LDStrength_vec = LD_Strength;
	ACD_RotStrength_vec = RM_Face;
	wM0Strength_vec = wM0_vec;
}

#endif // ABSORPTIONCD_HPP_INCLUDED
