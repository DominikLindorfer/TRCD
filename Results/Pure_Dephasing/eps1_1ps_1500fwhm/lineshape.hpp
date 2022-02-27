/*
 * lineshape.hpp
 *
 *  Created on: Jul 6, 2018
 *      Author: lindorfer
 */

#ifndef LINESHAPE_HPP_
#define LINESHAPE_HPP_

#include <string>
#include <cmath>
#include "Utilities.hpp"
//#include "correlation.hpp"
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <fftw3.h>
#include <complex>
#include <vector>
#include <cstdlib>
#include "four1.hpp"

using namespace std;

double J_WZ(double w){
	//-----Empirical Function of Spectral Density J(w) (2.16) in Wellenzahlen-----
    double huang = 1.0/1.3; //0.8/1.3; //1.0;

    double gamma = 53.04; //100fs

    if( w >= 0){
        return (huang*( 0.8/(5040.0 * 2 * pow(0.069e-3 / 6.582e-16 / 1.88496e11,4)) * pow(w,3) * exp(- sqrt(w / (0.069e-3 / 6.582e-16 / 1.88496e11))) + 0.5/(5040.0 * 2 * pow(0.24e-3 / 6.582e-16 / 1.88496e11,4)) * pow(w,3) * exp(- sqrt(w / (0.24e-3 / 6.582e-16 / 1.88496e11))) ));
    }

    if( w < 0){
        return 0;
    }
    else{
    	return -99999999;
    }

}

double J_WZ_ohm(double& w){

	double S_huang = 0.22;
	double J = 0;
	double lambda_c = 240;
	double gamma = 600;

	J = S_huang * lambda_c * gamma / (w * w + gamma * gamma);

	return J;
}

void Check_SpectralDensity(string Filename){

	ofstream out;
	out.open(Filename);

	for(double w = 0; w < 1000; w = w + 0.1){

		out << w << " " << J_WZ(w) << " " << J_WZ_ohm(w) << endl;
	}

	out.close();
}

double n_WZ(double w, double T){
	//-----Bose-Einstein Distribution (2.10) in Wellenzahlen-----
    if (w > 1e-6){
        return ( 1.0 / ( exp( (6.582 * pow(10.0,-16) * 1.88496 * pow(10.0,11) * w) / (8.617 * pow(10,-5) * T ) ) - 1.0 ));
    }
    else{
        return 1.0;
    }
}

double C_Re_WZ (double w, double T){
	//-----C_Re (4.15) in w [cm^-1]-----
//	return( 3.1415926535897932384626434 * w * w * J_WZ(w));
	return( 3.1415926535897932384626434 * w * w * ((1.0 + n_WZ(w,T)) * J_WZ(w) + n_WZ(-w,T) * J_WZ(-w)) );
}

double y_MNKL (int M, int N, int K, int L, Matrix<double>& evecs,  Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){
    //----- yMNKL from Eq. (4.9) -----
    //----- cm(M) etc. are the eigencoefficients, Rmn is the distance between the chlorophylls in question -----
    double sum = 0;
    double Rmn;
    Vektor<double> distance;

    for(int i = 0; i < (int)evecs.size1(); i++){
        for(int j = 0; j < (int)evecs.size2(); j++){
            distance = (dipolemat(1)(i) - dipolemat(1)(j));     //diolemat(1)(i) contains the centers-coordinates of the chlorophylls

            double norm = 0;
            norm = distance(0)*distance(0) + distance(1)*distance(1) + distance(2)*distance(2);

            if(norm < 1e-8){
                distance(0) = 0;
                distance(1) = 0;
                distance(2) = 0;
            }

//            Rmn = 1e-10 * norm_2 (distance);
            Rmn = 1e-10 * norm_2 (distance);
            sum += exp(- Rmn / parameter_vec(1)) * evecs(M,i) * evecs(N,i) * evecs (K,j) * evecs(L,j);

//            if(Rmn > 1e-5){
//            cout<<"Rmn= "<<Rmn<<endl;
//            cout<<"i = "<<i<<" j = "<<j<<" "<< /*exp(- Rmn / parameter_vec(1)) **/ evecs(M,i) * evecs(N,i) * evecs (K,j) * evecs(L,j)<<endl;
//            cout<<setprecision(15)<<evecs(M,i) * evecs(N,i) * evecs (K,j) * evecs(L,j)<<",";
//            }
        }
    }

    //errorchecking
    /*for(int i = 0; i < evecs.size1(); i++){
        cout<<evecs(L,i)<<",";
    }*/
    //cout<<evecs(0,0) * evecs(1,0) * evecs (2,0) * evecs(3,0)<<endl;
    //cout<<"Chl1: "<<dipolemat(1)(1)<<endl;
    //cout<<"Chl1: "<<dipolemat(1)(2)<<endl;
    //cout<<"InnerProd: "<<inner_prod(dipolemat(1)(1),dipolemat(1)(2))<<endl;
    //cout<<"Length: "<<dipolemat(1)(1)-dipolemat(1)(2)<<endl;
    //cout<<"Norm: "<<norm_2(dipolemat(1)(1)-dipolemat(1)(2))<<endl;
    return sum;
}

double y_MN (int M, int N, Matrix<double>& evecs){

	double sum = 0;

    for(int i = 0; i < (int)evecs.size1(); i++){

            sum += evecs(M,i) * evecs(M,i) * evecs(N,i) * evecs(N,i);
    }

    return sum;
}

double y_MN_Huang(int& M, int& N, Matrix<double>& evecs, Vektor<double>& Huang_Domain){

	double sum = 0;

    for(int i = 0; i < (int)evecs.size1(); i++){

            sum += evecs(M,i) * evecs(M,i) * evecs(N,i) * evecs(N,i) * Huang_Domain(i);
    }

    return sum;
}

double tau_WZ (int M, Vektor<double>& evals, Matrix<double>& evecs,  Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	double sum = 0;

	for(int K = 0; K < (int)evals.size(); K++){
        if( K != M){
            sum += y_MN(M, K, evecs) * C_Re_WZ(evals(M) - evals(K) , parameter_vec(0) );
        }
    }

    //double puredeph = 1.0/18.85;   //100000fs
    double puredeph = 1.0/0.1885;     //1000fs ----- 1000 10^-15 *(3 10^10)*2 Pi = 0.1885 in cm^-1
//    puredeph = 1.0/0.00231947;	//15fs
//    puredeph = 1.0/0.0188496;	//100fs
//    puredeph = 1.0/0.00942478;	//50fs
    return (sum + puredeph);
}

double tau_WZ_Huang (int M, Vektor<double>& evals, Matrix<double>& evecs,  Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& Huang_Domain){

	double sum = 0;

    for(int K = 0; K < (int)evals.size(); K++){
        if( K != M){

			sum += y_MN_Huang(M, K, evecs, Huang_Domain) * C_Re_WZ(evals(M) - evals(K), parameter_vec(0));
        }
    }

    //double puredeph = 1.0/18.85;   //100000fs
    double puredeph = 1.0/(0.3*0.1885);     //1000fs ----- 1000 10^-15 *(3 10^10)*2 Pi = 0.1885 in cm^-1
//    double puredeph = 0.0;

    return (sum + puredeph);
}

double wM0_WZ ( int M, Vektor<double>& evals, Matrix<double>& evecs, Vektor<double>& parameter_vec, Vektor<double>& CIm, Vektor< Vektor< Vektor<double> > >& dipolemat){

    double sum = 0;
    double wM0;

    double te = 8192.0/512.0;
    double ntime = 8192.0*32;
    double twopi = 3.141592 * 2.0;
    int suche;
    double CImags;

    for( int K = 0; K < (int)evals.size(); K++){
        if( K != M){

            suche = (evals(M) - evals(K) + (twopi/te*(ntime/2.0 -1)))*te/(twopi);
            CImags = CIm(suche);

            sum += y_MN(M, K, evecs) * CImags;
        }
    }

    wM0 = evals(M) - y_MN(M, M, evecs) * 101.673 / 1.3 + sum;

    return( wM0 );
}

double wM0_WZ_Huang ( int M, Vektor<double>& evals, Matrix<double>& evecs, Vektor<double>& parameter_vec, Vektor<double>& CIm, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double>& Huang_Domain){

    double sum = 0;
    double wM0;

    double te = 8192.0/512.0;
    double ntime = 8192.0*32;
    double twopi = 3.141592 * 2.0;
    int suche;
    double CImags;

    for( int K = 0; K < (int)evals.size(); K++){
        if( K != M){

			suche = (evals(M) - evals(K) + (twopi / te * (ntime / 2.0 - 1))) * te / (twopi);
            CImags = CIm(suche);

			sum += y_MN_Huang(M, K, evecs, Huang_Domain) * CImags;
        }
    }

//    cout << "sum: " << sum << endl;
//    cout << "E_lambda: " << y_MN_Huang(M, M, evecs, Huang_Domain) * 102 << endl;

	wM0 = evals(M) - y_MN_Huang(M, M, evecs, Huang_Domain) * 101.673/1.3 + sum;
    //wM0 = evals(M) - y_MNKL(M,M,M,M,evecs,parameter_vec,dipolemat) * 101.673;         //Julian

    return( wM0 );
}

double y_2N2K (int N2, int K2, Matrix<double>& evec2N, Matrix<double>& Exc2_Map, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){
	//-----Eq. 4.55 in the R/M Paper-----

    //-----Exc2_Map maps from 1-Exciton states to the corresponding 2-Exciton state-----
    //-----Exc2_Map.size() is the 1-Exciton Manifold size. Mapping is done when calculating 'sum'-----

	double sum = 0;

    for(int m = 0 ; m < (int)Exc2_Map.size1(); m++){
    	for(int n = 0; n < (int)Exc2_Map.size1(); n++){

    		for(int k = 0; k < (int)Exc2_Map.size1(); k++){
    			for(int l = 0; l < (int)Exc2_Map.size1(); l++){

    				if(m > n){
    					if(k > l){

    						if(m == k){
    							sum += evec2N(N2, Exc2_Map(m,n)) * evec2N(N2, Exc2_Map(m,n)) * evec2N(K2, Exc2_Map(k,l)) * evec2N(K2, Exc2_Map(k,l));
    						}
    						if(m == l){
    							sum += evec2N(N2, Exc2_Map(m,n)) * evec2N(N2, Exc2_Map(m,n)) * evec2N(K2, Exc2_Map(k,l)) * evec2N(K2, Exc2_Map(k,l));
    						}
    						if(n == k){
    							sum += evec2N(N2, Exc2_Map(m,n)) * evec2N(N2, Exc2_Map(m,n)) * evec2N(K2, Exc2_Map(k,l)) * evec2N(K2, Exc2_Map(k,l));
    						}
    						if(n == l){
    							sum += evec2N(N2, Exc2_Map(m,n)) * evec2N(N2, Exc2_Map(m,n)) * evec2N(K2, Exc2_Map(k,l)) * evec2N(K2, Exc2_Map(k,l));
    						}

    					}
    				}

    			}
    		}
    	}
    }

    return sum;
}

double y_2N2KMM (int N2, int M, Matrix<double>& evec2N, Matrix<double>& evecs, Matrix<double>& Exc2_Map, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){
	//-----Eq. 4.55 in the R/M Paper-----
    double sum = 0;

    for(int m = 0 ; m < (int)Exc2_Map.size1(); m++){
    	for(int n = 0; n < (int)Exc2_Map.size1(); n++){

    		for (int k = 0; k < (int) evecs.size1(); k++) {

				if (m > n) {

					if(m == k){
						sum += evec2N(N2, Exc2_Map(m, n))* evec2N(N2, Exc2_Map(m, n)) * evecs(M, k) * evecs(M, k);
					}

					if(n == k){
						sum += evec2N(N2, Exc2_Map(m, n))* evec2N(N2, Exc2_Map(m, n)) * evecs(M, k) * evecs(M, k);
					}
				}
    		}
    	}
    }

    return sum;
}

double w2N0_WZ ( int N2, Vektor<double>& evals2N, Matrix<double>& evecs2N, Matrix<double>& Exc2_Map, Vektor<double>& parameter_vec, Vektor<double>& CIm, Vektor< Vektor< Vektor<double> > >& dipolemat){

    double w2N0 = evals2N(N2) - y_2N2K(N2, N2, evecs2N, Exc2_Map, parameter_vec,dipolemat) * 101.673 / 1.3;
    return( w2N0 );
}

double w2NM_WZ (double w2N0, double wM0){

    return( w2N0 - wM0 );
}

double tau2N_WZ(int N2, Vektor<double>& evals2N, Matrix<double>& evecs2N, Matrix<double>& Exc2_Map, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	double sum = 0;

    for(int K2 = 0; K2 < (int)evals2N.size(); K2++){
        if( K2 != N2){

            sum += y_2N2K(N2, K2, evecs2N, Exc2_Map, parameter_vec, dipolemat) * C_Re_WZ(evals2N(N2) - evals2N(K2) , parameter_vec(0) );
        }
    }

    double puredeph = 1.0/0.1885;     //1000fs
    puredeph = 1.0/0.00231947;
//    return sum;
    return (sum + puredeph);
}

double tau2NM_WZ(double tau2N, double tauM){

	return(tau2N + tauM);
}

void FT_Gt_WZ(Vektor<double>& GReals, Vektor<double>& GImags, double te, double ntime, Vektor<double> parameter_vec){  //mit four1-Routine

    double N = ntime;
    double delta = te / (ntime -1);
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;

    Vektor<double> spec(2 * (int)ntime + 1);

    //Data alignment in in[i][j] accoriding to four1 from above
    for(int i = 1; i <= ntime/2 -1; i++){
         om = twopi * i / delta / ntime;                        //wir sind hier im k-Raum; deswegen muss spacing so gewählt werden

         spe = (1.0 + n_WZ(om,parameter_vec(0))) * J_WZ(om);
         spe1 = n_WZ(om,parameter_vec(0)) * J_WZ(om);

         spec(2*i + 1) = spe;
         spec(2*(int)ntime - 2*i+1) = spe1;
         spec(2*i+2) = 0.;
         spec(2*(int)ntime - 2*i+2) = 0.;
    }

    spec(1) = 0.0;
    spec(2) = 0.0;
    spec((int)ntime +1) = 0;
    spec((int)ntime +2) = 0;

    int is = -1;
    double corr_real, corr_imag;

    //cout<<"G(t) vor FFT"<<endl;

    four1_FFTW(spec, ntime, is);
    //Output for Errorcorrections
//    ofstream myfile1;
//    myfile1.open ("test_G(t).dat");

    GReals.resize(ntime);
    GImags.resize(ntime);
    //cout<<"G(t) nach FFT"<<endl;
    Vektor<double> corr(2 * (int)ntime + 1);

    for(int i = 1; i <= ntime; i++){

        corr_real = spec(2*i-1) / ntime * twopi / delta;
        corr_imag = spec(2*i) / ntime * twopi / delta;

        GReals(i-1) = corr_real;
        GImags(i-1) = corr_imag;
//        myfile1 <<setprecision(15)<< (i-1)*delta /*/(pow(10,-15))*/ <<" "<<GReals(i-1)<<" "<<GImags(i-1)<<endl;    //<<" "<<exp(-i*delta)
    }
    //cout<<"GReals(0) = "<<GReals(0)<<endl;
//    myfile1.close();
//    cout<<"G(t) am Ende von FFT"<<endl;
}

void DM_WZ(int M, Vektor<double>& DMReals,Vektor<double>& GReal, Vektor<double>& GImag, double te, double ntime, Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, fftw_complex* in, fftw_complex* out){

	//-----Initialize Parameter-----
	Vektor<double> GReals;
    Vektor<double> GImags;

    GReals = GReal * y_MN(M, M, evecs);
    GImags = GImag * y_MN(M, M, evecs);

    double GM0 = GReals(0);
    double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
    double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);

    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //-----Fill Up the FFTs-----
    for(int i = 0; i < ntime/2 ; i++){

         in[i][0] = exp(GReals(i) - (i+1) * delta * TauM - GM0) * cos(GImags(i));
         in[i][1] = exp(GReals(i) - (i+1) * delta * TauM - GM0) * sin(GImags(i));

    }

    for(int i = ntime-1; i >= ntime/2; i = i-1){

        in[i][0] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * cos(GImags(i));
        in[i][1] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * sin(GImags(i));
    }

    fftw_plan p;
    #pragma omp critical (make_plan)
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    fftw_execute_dft(p,in,out);
    fftw_destroy_plan(p);

    DMReals.resize(ntime);

    double spectral;

    spectral = in[0][0] * delta;

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

    //-----Cut-out proper w-Range -----
    //----- ntime/2 is the w=0 point in DMReals -----
    //----- te/twopi * wM0 shifts the curve to it's actual position after the FFT-----
    int upperlim = ntime/2.0  + te/twopi * (1.0/(300*1e-7) - wM0);
    int lowerlim = ntime/2.0  + te/twopi * (1.0/(800*1e-7) - wM0);

//    cout << "upperlim : "<<upperlim<<endl;
//    cout << "lowerlim : "<<lowerlim<<endl;
    //cout<<"Difference Upper/Lower: "<<upperlim-lowerlim<<endl;          //Check here for errors, (int) lower/upperlim might not be equal for all States M!
    //exit(1);

    Vektor<double> Temp(upperlim - lowerlim);

    int tempsize = (int)Temp.size();
    for(int i = 0; i < (upperlim-lowerlim); i++){

                Temp(tempsize-1 - i) = DMReals(lowerlim + i);
    }

    DMReals.resize(13262);

    for(int i = 0; i < 13262; i++){
            DMReals(i) = Temp(i);
    }
}

void DM_Emission_WZ(int M, Vektor<double>& DMReals,Vektor<double>& GReal, Vektor<double>& GImag, double te, double ntime, Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, fftw_complex* in, fftw_complex* out){

	//Initialize Parameter
    Vektor<double> GReals;
    Vektor<double> GImags;

    GReals = GReal * y_MN(M, M, evecs);
    GImags = GImag * y_MN(M, M, evecs);

    double GM0 = GReals(0);
    double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
    double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);

    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //-----Fill Up the FFTs-----
    for(int i = 0; i < ntime/2 ; i++){

         in[i][0] = exp( (GReals(i) - GM0) - (i+1) * delta * TauM ) * cos(GImags(i));
         in[i][1] = exp( (GReals(i) - GM0) - (i+1) * delta * TauM ) * sin(GImags(i));

    }

    for(int i = ntime-1; i >= ntime/2; i = i-1){

        in[i][0] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * cos(GImags(i));
        in[i][1] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * sin(GImags(i));
    }

    fftw_plan p;
    #pragma omp critical (make_plan)
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute_dft(p,in,out);
    fftw_destroy_plan(p);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

    //-----Cut-out proper w-Range -----
    //----- ntime/2 is the w=0 point in DMReals -----
    //----- te/twopi * wM0 shifts the curve to it's actual position after the FFT-----

    int upperlim = ntime/2.0  + te/twopi * (1.0/(300*1e-7) - wM0);      //approximatively this  works [600,730]!!!
    int lowerlim = ntime/2.0  + te/twopi * (1.0/(800*1e-7) - wM0);

    Vektor<double> Temp(upperlim - lowerlim);
    int tempsize = (int)Temp.size();

    for(int i = 0; i < (upperlim-lowerlim); i++){

                Temp(tempsize-1 - i) = DMReals(lowerlim + i);
    }

    DMReals.resize(13262);

    for(int i = 0; i < 13262; i++){
            DMReals(i) = Temp(i);
    }
}

void DM_Emission_WZ_Huang(int M_D1, int dom_a, Vektor<double>& DMReals,Vektor<double>& GReal, Vektor<double>& GImag, double te, double ntime, Vektor<double>& CIm, Vektor< Matrix<double> >& evecs_vec, Vektor< Vektor<double> >& evals_vec, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor< Vektor<double> >& Huang_Domains, fftw_complex* in, fftw_complex* out, fftw_plan& p){

	//Initialize Parameter
    Vektor<double> GReals;
    Vektor<double> GImags;

    GReals = GReal * y_MN_Huang(M_D1, M_D1, evecs_vec(dom_a), Huang_Domains(dom_a));
    GImags = GImag * y_MN_Huang(M_D1, M_D1, evecs_vec(dom_a), Huang_Domains(dom_a));

    double GM0 = GReals(0);
	double TauM = tau_WZ_Huang(M_D1, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, dipolemat, Huang_Domains(dom_a));
	double wM0 = wM0_WZ_Huang(M_D1, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, CIm, dipolemat, Huang_Domains(dom_a));

    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //-----Fill Up the FFTs-----
    for(int i = 0; i < ntime/2 ; i++){

         in[i][0] = exp( (GReals(i) - GM0) - (i+1) * delta * TauM ) * cos(GImags(i));
         in[i][1] = exp( (GReals(i) - GM0) - (i+1) * delta * TauM ) * sin(GImags(i));

    }

    for(int i = ntime-1; i >= ntime/2; i = i-1){

        in[i][0] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * cos(GImags(i));
        in[i][1] = exp(GReals(i) + (i+1-ntime) * delta * TauM - GM0) * sin(GImags(i));
    }

    #pragma omp critical (make_plan)
    {
//        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute_dft(p,in,out);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

    //-----Cut-out proper w-Range -----
    //----- ntime/2 is the w=0 point in DMReals -----
    //----- te/twopi * wM0 shifts the curve to it's actual position after the FFT-----

    int upperlim = ntime/2.0  + te/twopi * (1.0/(300*1e-7) - wM0);      //approximatively this  works [600,730]!!!
    int lowerlim = ntime/2.0  + te/twopi * (1.0/(800*1e-7) - wM0);

    Vektor<double> Temp(upperlim - lowerlim);

    int tempsize = (int)Temp.size();

    for(int i = 0; i < (upperlim-lowerlim); i++){

                Temp(tempsize-1 - i) = DMReals(lowerlim + i);
    }

    DMReals.resize(13262);    //Dirty Solution to get a hold of Error above, cutting out 1-2 values

    for(int i = 0; i < 13262; i++){
            DMReals(i) = Temp(i);
    }
}

void DM_ESA_WZ(int N2, int M, Vektor<double>& DMReals,Vektor<double> GReal, Vektor<double> GImag, double te, double ntime, Vektor<double> CIm, Matrix<double> evecs, Vektor<double> evals, Vektor<double> evals2N, Matrix<double> evecs2N, Matrix<double> Exc2_Map, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat, fftw_complex* in, fftw_complex* out){

	//Initialize Parameter
    Vektor<double> GReals;
    Vektor<double> GImags;

    GReals = GReal * ( y_2N2K(N2, N2, evecs2N, Exc2_Map, parameter_vec, dipolemat) + y_MNKL(M,M,M,M,evecs,parameter_vec,dipolemat) - 2 * y_2N2KMM(N2, M, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) );
    GImags = GImag * ( y_2N2K(N2, N2, evecs2N, Exc2_Map, parameter_vec, dipolemat) + y_MNKL(M,M,M,M,evecs,parameter_vec,dipolemat) - 2 * y_2N2KMM(N2, M, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) );

    double GM0 = GReals(0);

    double TauM  = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
    double Tau2N = tau2N_WZ(N2, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat);
    double Tau2NM = Tau2N + TauM;

	double wM0  = wM0_WZ(M, evals, evecs, parameter_vec, CIm, dipolemat);
	double w2N0 = w2N0_WZ(N2, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat);
	double w2NM = w2N0 - wM0;

    //int N = ntime;
    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //*****Fill Up the FFTs*****
	for (int i = 0; i < ntime / 2; i++) {

		in[i][0] = exp((GReals(i) - GM0) - (i + 1) * delta * Tau2NM) * cos(GImags(i));
		in[i][1] = exp((GReals(i) - GM0) - (i + 1) * delta * Tau2NM) * sin(GImags(i));

	}

	for (int i = ntime - 1; i >= ntime / 2; i = i - 1) {

		in[i][0] = exp((GReals(i) - GM0) + (i + 1 - ntime) * delta * Tau2NM) * cos(GImags(i));
		in[i][1] = exp((GReals(i) - GM0) + (i + 1 - ntime) * delta * Tau2NM) * sin(GImags(i));
	}

    fftw_plan p;
    #pragma omp critical (make_plan)
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    fftw_execute_dft(p,in,out);
    fftw_destroy_plan(p);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
        //myfile1 <<setprecision(15)<< i*delta <<" "<<GReals(i)<<" "<<GImags(i)<<endl;    //<<" "<<exp(-i*delta)
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

//    ofstream myfile;
//    myfile.open("test_Dm_forward.dat");
//    for(int i = 0; i<ntime; i++){
//
//        myfile<<(i-ntime/2.0)*twopi/te<<" "<<DMReals(i)<<endl;
//    }
//    myfile.close();

    //*****Cut-out proper w-Range (this is a dirty solution!)*****
    //***** ntime/2 is the w=0 point in DMReals *****
    //***** te/twopi * wM0 shifts the curve to it's actual position after the FFT*****

    int upperlim = ntime/2.0  + te/twopi * (1.0/(300*1e-7) - w2NM);      //approximatively this  works [600,730]!!!
    int lowerlim = ntime/2.0  + te/twopi * (1.0/(800*1e-7) - w2NM);

//    //-----DEBUG-----
//	cout <<"DEBUG --" << "wM0: " << wM0 << " w2N0: " << w2N0 << " w2NM: " << w2NM << " upperlim : " << upperlim << "   lowerlim : " << lowerlim << endl;

    //cout<<"Difference Upper/Lower: "<<upperlim-lowerlim<<endl;          //Check here for errors, (int) lower/upperlim might not be equal for all States M!
    //241857
    //exit(1);
    //cout<<"DMReals Lower , Upper: "<<DMReals(lowerlim)<<" "<<DMReals(upperlim)<<endl;

    //Cut-out proper w-Range (this is a dirty solution!)
    Vektor<double> Temp(upperlim - lowerlim);
//    cout<<"Size of Temp: "<<Temp.size()<<endl;
//    cout<<"ExcitonState M : "<<M<<endl;
//    cout<<"Eval M : "<<evals(M)<<endl;
//    cout<<"wM0 : "<<wM0<<endl;

    int tempsize = (int)Temp.size();

    for(int i = 0; i < (upperlim-lowerlim); i++){

                Temp(tempsize-1 - i) = DMReals(lowerlim + i);
                //myfile<<((i+lowerlim)-ntime/2.0)*twopi/te<<" "<<Temp(i)<<endl;
                //Temp(Temp.size()-1 - i) = DMReals[lowerlim+i];    //works for [600,730]
                //myfile<<((i+lowerlim)-ntime/2.0)*twopi/te<<" "<<Temp(i)<<endl;
    }

//    cout << "Before Resize:" <<endl;

    //myfile.close();4194304
    //1697652   53052
    DMReals.resize(13262);    //Dirty Solution to get a hold of Error above, cutting out 1-2 values

    for(int i = 0; i < 13262; i++){   //original: 1697652  //for 8192.0/256.0 ,8192.0 * 32.0 :  106102
            DMReals(i) = Temp(i);
    }

//    cout << "After Resize" <<endl;

    //myfile1.close();

//
//    fftw_destroy_plan(p);
//    fftw_free(in);
//    fftw_free(out);
}

void DM_WZ_Markov(int M, Vektor<double>& DMReals,Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	//-----Initialize Parameter-----
    double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat); // 2 * M_PI *
    double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);

    long upperlim = 1.0/(400*1e-7);
    long lowerlim = 1.0/(1000*1e-7);

    DMReals.resize(upperlim - lowerlim);

    for(long w = lowerlim; w < upperlim; w++){

    	DMReals(w - lowerlim) = TauM / ((w - wM0)*(w - wM0) + TauM * TauM);
    }

}

void DM_Emission_WZ_Markov(int M, Vektor<double>& DMReals,Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	//-----Initialize Parameter-----
    double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat); // 2 * M_PI *
    double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);

    long upperlim = 1.0/(400*1e-7);
    long lowerlim = 1.0/(1000*1e-7);

    DMReals.resize(upperlim - lowerlim);

    for(long w = lowerlim; w < upperlim; w++){

    	DMReals(w - lowerlim) = TauM / ((-w + wM0)*(-w + wM0) + TauM * TauM);
    }

}

void DM_ESA_WZ_Markov(int N2, int M, Vektor<double>& DMReals, Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& evals2N, Matrix<double>& evecs2N, Matrix<double> Exc2_Map, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	//-----Initialize Parameter-----
    double TauM  = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
    double Tau2N = tau2N_WZ(N2, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat);
    double Tau2NM = Tau2N + TauM;

	double wM0  = wM0_WZ(M, evals, evecs, parameter_vec, CIm, dipolemat);
	double w2N0 = w2N0_WZ(N2, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat);
//	double w2NM = w2N0 - wM0 + 90;
	double w2NM = w2N0 - wM0;

	cout << "N2: " << N2 << "  M:" << M << " wM0: " << wM0 << " w2NM: " << w2NM << " w2N0: " << w2N0 << " TauM: " << TauM << " Tau2NM: " << Tau2NM << " y2N2K: " << y_2N2K(N2, N2, evecs2N, Exc2_Map, parameter_vec,dipolemat) << endl;

    long upperlim = 1.0/(400*1e-7);
    long lowerlim = 1.0/(1000*1e-7);

    DMReals.resize(upperlim - lowerlim);

    for(long w = lowerlim; w < upperlim; w++){

    	DMReals(w - lowerlim) = Tau2NM / ((w - w2NM)*(w - w2NM) + Tau2NM * Tau2NM);
    }
}

void DM_WZ_Markov_old(int M, Vektor<double>& DMReals, Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat){
    //Initialize Parameter

	double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
	double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);

	//-----lower_lim is in lambda -> convert to w values-----
	DMReals.resize(13262);
	long u = DMReals.size();
	long lower_lim = 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0)));

	for(long i = 0; i < DMReals.size(); i++){

		double w_lambda = 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - i * (2 * 3.1415926) / (4.0)));
		double w_WZ = 1.0 / (w_lambda * 1e-7);

		DMReals(i) = TauM / (pow((w_WZ - wM0), 2) + pow(TauM, 2));
	}
}

double PM0_WZ(int M, Vektor<double>& DMReals,Vektor<double>& GReal, Vektor<double>& GImag, double te, double ntime, Vektor<double>& CIm, Matrix<double>& evecs, Vektor<double>& evals, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, fftw_complex* in, fftw_complex* out, double pulse_center_omega, double TauP){

	//Initialize Parameter
    Vektor<double> GReals;
    Vektor<double> GImags;

    GReals = GReal * y_MN(M,M,evecs);
    GImags = GImag * y_MN(M,M,evecs);

    double GM0 = GReals(0);

    double TauM  = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);
	double wM0  = wM0_WZ(M, evals, evecs, parameter_vec, CIm, dipolemat);

    //int N = ntime;
    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //-----Fill Up the FFTs-----
	for (int i = 0; i < ntime / 2; i++) {

		in[i][0] = exp((GReals(i) - GM0) - ((i + 0) * delta) * ((i + 0) * delta) / (4 * TauP * TauP)) * cos(GImags(i));
		in[i][1] = exp((GReals(i) - GM0) - ((i + 0) * delta) * ((i + 0) * delta) / (4 * TauP * TauP)) * sin(GImags(i));

	}

	for (int i = ntime - 1; i >= ntime / 2; i = i - 1) {

		in[i][0] = exp((GReals(i) - GM0) - ((i + 0 - ntime) * delta) * ((i + 0 - ntime) * delta) / (4 * TauP * TauP)) * cos(GImags(i));
		in[i][1] = exp((GReals(i) - GM0) - ((i + 0 - ntime) * delta) * ((i + 0 - ntime) * delta) / (4 * TauP * TauP)) * sin(GImags(i));
	}

	fftw_plan p;
    #pragma omp critical (make_plan)
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute_dft(p,in,out);
    fftw_destroy_plan(p);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
        //myfile1 <<setprecision(15)<< i*delta <<" "<<GReals(i)<<" "<<GImags(i)<<endl;    //<<" "<<exp(-i*delta)
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

////  -----Output DMReals-----
//	ofstream myfile;
//	myfile.open("test_P0_GB.dat");
//	for (int i = 0; i < ntime; i++) {
//
//		myfile << (i - ntime / 2.0) * twopi / te << " " << DMReals(i) << endl;
//	}
//	myfile.close();

//    cout << "wM0 = " << wM0 << endl;

    //-----Give back proper Integral-value at frequency (omega - wM0)-----
	int excitation_value = ntime / 2.0 + te / twopi * (1.0 / (pulse_center_omega * 1e-7) - wM0);
	return ( 1e4 * DMReals(excitation_value) );

}

double GF_kMK_WZ(int M_D1, int K_D2, int dom_a, int dom_b, Vektor<double>& DMReals, Vektor<double>& GReal, Vektor<double>& GImag, double te, double ntime, Vektor<double>& CIm, Vektor< Matrix<double> >& evecs_vec, Vektor< Vektor<double> >& evals_vec, Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor< Vektor<double> >& Huang_Domains, fftw_complex* in, fftw_complex* out, fftw_plan& p){

	//Initialize Parameter
    Vektor<double> GRealsM;
    Vektor<double> GImagsM;

    Vektor<double> GRealsK;
    Vektor<double> GImagsK;

    GRealsM = GReal * y_MN_Huang(M_D1, M_D1, evecs_vec(dom_a), Huang_Domains(dom_a));
	GImagsM = GImag * y_MN_Huang(M_D1, M_D1, evecs_vec(dom_a), Huang_Domains(dom_a));

	GRealsK = GReal * y_MN_Huang(K_D2, K_D2, evecs_vec(dom_b), Huang_Domains(dom_b));
	GImagsK = GImag * y_MN_Huang(K_D2, K_D2, evecs_vec(dom_b), Huang_Domains(dom_b));

    double GM0 = GRealsM(0);
    double GK0 = GRealsK(0);

    double wM0 = wM0_WZ_Huang(M_D1, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, CIm, dipolemat, Huang_Domains(dom_a));
    double wK0 = wM0_WZ_Huang(K_D2, evals_vec(dom_b), evecs_vec(dom_b), parameter_vec, CIm, dipolemat, Huang_Domains(dom_b));

//    if(wM0 < wK0){
//    	return 0;
//    }

	double TauM = tau_WZ_Huang(M_D1, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, dipolemat, Huang_Domains(dom_a));
	double TauK = tau_WZ_Huang(K_D2, evals_vec(dom_b), evecs_vec(dom_b), parameter_vec, dipolemat, Huang_Domains(dom_b));

//	cout << "GF Function: " << endl;
//	cout << M_D1 << " " << K_D2 << " " << dom_a << " " << dom_b << endl;
//	cout << wM0 << " " << wK0 << endl;
//	cout << TauM << " " << TauK << endl;
//	cout << GM0 << " " << GK0 << endl;
//	cout << y_MN_Huang(M_D1, M_D1, evecs_vec(dom_a), Huang_Domains(dom_a)) << endl;

    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //-----Fill Up the FFTs-----

	for (int i = 0; i < ntime / 2; i++) {

		in[i][0] = exp((GRealsM(i)- GM0) + (GRealsK(i)- GK0) - (i + 1) * delta * (TauM + TauK) ) * cos(GImagsM(i) + GImagsK(i));
		in[i][1] = exp((GRealsM(i)- GM0) + (GRealsK(i)- GK0) - (i + 1) * delta * (TauM + TauK) ) * sin(GImagsM(i) + GImagsK(i));
	}

	for (int i = ntime - 1; i >= ntime / 2; i = i - 1) {

		in[i][0] = exp((GRealsM(i)- GM0) + (GRealsK(i)- GK0) + (i + 1 - ntime) * delta * (TauM + TauK)) * cos(GImagsM(i) + GImagsK(i));
		in[i][1] = exp((GRealsM(i)- GM0) + (GRealsK(i)- GK0) + (i + 1 - ntime) * delta * (TauM + TauK)) * sin(GImagsM(i) + GImagsK(i));
	}

    #pragma omp critical (make_plan)
    {
        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
//        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute_dft(p,in,out);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){
        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
    }

    for(int i = ntime/2.0; i< ntime; i++){
        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }

////  -----Output DMReals-----
//	ofstream myfile;
//	myfile.open("test_GeneralizedFoerster.dat");
//	cout << "Print Generalized Foerster to File!" << endl;
//	for (int i = 0; i < ntime; i++) {
//
//		myfile << (i - ntime / 2.0) * twopi / te << " " << DMReals(i) << endl;
//	}
//	myfile.close();
//    cout << "wM0 = " << wM0 << endl;
//    cout << "wK0 = " << wK0 << endl;

    //-----Give back proper Integral-value at frequency (wK0 - wM0)-----
	int excitation_value = ntime / 2.0 + te / twopi * (wM0 - wK0);

//	cout << DMReals(excitation_value) << endl;
	return ( DMReals(excitation_value) );

}

void DM_Vorfaktor(Vektor<double>& DMReals, double te, double ntime, fftw_complex* in, fftw_complex* out, fftw_plan& p){

	//Initialize Parameter
    //int N = ntime;
    double N = ntime;    //gives the size of the Vectors you put into FFTW and get out of it; The bigger this is the more accurate the results
    double delta = te / (ntime -1);             //this sets the interval for frequency / time
    double twopi = 2.0 * 3.1415926535897932384626434;
    double spe,spe1,om;
    double time;
    double re,im;

    //*****Fill Up the FFTs*****
    for(int i = 0; i < ntime/2 ; i++){

         in[i][0] = exp(-(i * delta) * (i * delta));
         in[i][1] = 0.0;

    }

    for(int i = ntime-1; i >= ntime/2; i = i-1){

        in[i][0] = exp( -((i-ntime) * delta) * ((i-ntime) * delta) );
        in[i][1] = 0.0;
    }

    #pragma omp critical (make_plan)
    {
//        p = fftw_plan_dft_1d(N, in, out, FFTW_BACKWARD, FFTW_ESTIMATE);
        p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    fftw_execute_dft(p,in,out);

    DMReals.resize(ntime);

    for(int i = 0; i < ntime/2.0; i++){

        DMReals(ntime/2.0 + i) = out[i][0] * delta / twopi;
        //myfile1 <<setprecision(15)<< i*delta <<" "<<GReals(i)<<" "<<GImags(i)<<endl;    //<<" "<<exp(-i*delta)
    }

    for(int i = ntime/2.0; i< ntime; i++){

        DMReals(i - ntime/2.0) = out[i][0] * delta / twopi;
    }


    ofstream myfile;
	myfile.open("test_LorentzExp_Vorfaktor_FFT.dat");
	for (int i = 0; i < ntime; i++) {

		myfile << (i - ntime / 2.0) * twopi / te << " " << DMReals(i) << endl;
	}
	myfile.close();
}

double DM_GB(double w, int M, Vektor<double> CIm, Matrix<double> evecs, Vektor<double> evals, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat){

    double TauM = tau_WZ(M,evals,evecs,parameter_vec,dipolemat);		//Initialize Parameters
    double wM0 = wM0_WZ(M,evals,evecs,parameter_vec,CIm,dipolemat);		//TauM & wMO

//    cout << "TauM :" << TauM << endl << "wM0: " << wM0 << endl;

    return ( TauM / ((w - wM0)*(w - wM0) + TauM * TauM));
}

#endif /* LINESHAPE_HPP_ */
