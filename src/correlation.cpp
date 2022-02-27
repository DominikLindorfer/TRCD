#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include "four1.hpp"
//#include "lineshape.hpp"
#include "Utilities.hpp"

using namespace std;

fftw_plan p1_FFTW;
fftw_plan p2_FFTW;

/*void correlation(Vektor<double> parameter_vec, Vektor<double>& CIm){

    double te = 8192.0/512.0;   //8192.0 /256.0; //works
    double ntime = 8192.0*32;    //8192.0 * 32;   //works
    double delta = te / (ntime -1);
    double om,spe,spe1,spectral;
    double twopi = 3.141592 * 2.0;

    double spec[2 * (int)ntime + 1];

    for(int i = 1; i <= ntime / 2.0 -1; i++){
         om = twopi * i / delta / ntime;
         //spe = 1.0 / (1.0 + om*om) / (0.5 * twopi);
         //spe1 = spe;
         spe = om*om * (1.0 + n_WZ(om,parameter_vec(0))) * J_WZ(om);
         spe1 = om*om * n_WZ(om,parameter_vec(0)) * J_WZ(om);

         spec[2*i + 1] = spe;
         spec[2*(int)ntime - 2*i+1] = spe1;
         spec[2*i+2] = 0.;
         spec[2*(int)ntime - 2*i+2] = 0.;
    }

    spec[1] = 0.0;  //1.0 / (0.5*twopi);
    spec[2] = 0.0;
    spec[(int)ntime +1] = 0;
    spec[(int)ntime +2] = 0;

    int is = -1;
    double corr_real, corr_imag;

    four1_FFTW(spec, ntime, is);

    //ofstream myfile;
    //myfile.open ("corr.dat");


    double *corr = 0;
    corr = new double [2 * (int)ntime +1];

    for(int i = 1; i <= ntime; i++){

        corr_real = spec[2*i-1] / ntime * twopi / delta;
        corr_imag = spec[2*i] / ntime * twopi / delta;
        //myfile <<setprecision(15)<< (i-1)*delta/(twopi*3.0*pow(10,10)) <<" "<<corr_real<<" "<<corr_imag<<" "<<endl;
        corr[2*i-1] = corr_real;
        corr[2*i] = corr_imag;
    }
    //myfile.close();

    double re, im;
    for(int i = 1; i<=ntime/2.0 ; i++){
            re = corr[2*i-1];
            im = corr[2*i];
            spec[2*i-1] = re;
            spec[2*i] = im;
    }

    for(int i = ntime/2.0 +1; i<=ntime; i++){
            spec[2*i-1] = 0.0;
            spec[2*i] = 0.0;
    }

    is = 1;
    four1_FFTW(spec, ntime, is);

    ofstream cimout;
    cimout.open ("test_C_IM.dat");

    //Prepare for passing to &vectors
    int k = 0;  //will go from 0 to ntime-1
    //CRe.resize((int) (ntime));
    CIm.resize((int) (ntime));

    //w-Frame from [-ntime/2.0 -1 , ntime/2.0 -1]
    //for ntime = 8192*32 -> [-51471.5 , 51471.5]  which is maybe a bit too rough for current problems

    for(int i = ntime/2.0 -1; i>=1;i = i-1){
            cimout<<-twopi*i/te<<" "<<spec[2*(int)ntime-2*i+1]*delta<<" "<<spec[2*(int)ntime-2*i+2]*delta<<" "<<endl;
            //CRe(k) = spec[2*(int)ntime-2*i+1]*delta;
            CIm(k) = spec[2*(int)ntime-2*i+2]*delta;
            k++;
    }

    for(int i = 0; i<=ntime/2-1;i++){
            cimout<<twopi*i/te<<" "<<spec[2 * i +1] * delta<<" "<<spec[2 * i +2] *delta<<" "<<endl;
            //CRe(k) = spec[2 * i +1]*delta;
            CIm(k) = spec[2 * i +2]*delta;
            k++;
    }
    cimout.close();
    delete[] corr;
    corr = 0;
}*/

double J_WZ(double w);
double n_WZ(double w, double T);
double y_MNKL(int M, int N, int K, int L, Matrix<double>& evecs,  Vektor<double>& parameter_vec, Vektor< Vektor< Vektor<double> > >& dipolemat);

void correlation(Vektor<double> parameter_vec, Vektor<double>& CIm){  //Ublas-Vektoren

    double te = 8192.0/512.0;   //8192.0 /256.0; //works
    double ntime = 8192.0*32;    //8192.0 * 32;   //works
    double delta = te / (ntime -1);
    double om,spe,spe1,spectral;
    double twopi = 3.141592 * 2.0;

    Vektor<double> spec(2 * (int)ntime + 1);

    for(double i = 1; i <= ntime / 2.0 -1; i++){
         om = twopi * i / delta / ntime;

         spe = om*om * (1.0 + n_WZ(om,parameter_vec(0))) * J_WZ(om);
         spe1 = om*om * n_WZ(om,parameter_vec(0)) * J_WZ(om);

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

    four1_FFTW(spec, ntime, is);

    //ofstream myfile;
    //myfile.open ("test_C(t).dat");

    Vektor<double> corr(2 * (int)ntime + 1);

    for(double i = 1; i <= ntime; i++){

        corr_real = spec(2*i-1) / ntime * twopi / delta;
        corr_imag = spec(2*i) / ntime * twopi / delta;
        //myfile <<setprecision(15)<< (i-1)*delta<<" "<<corr_real<<" "<<corr_imag<<" "<<endl;
        corr(2*i-1) = corr_real;
        corr(2*i) = corr_imag;
    }
    //myfile.close();

    for(double i = 1; i<=ntime/2.0 ; i++){
            spec(2*i-1) = corr(2*i - 1);
            spec(2*i) = corr(2*i);
    }

    for(double i = ntime/2.0 +1; i<=ntime; i++){
            spec(2*i-1) = 0.0;
            spec(2*i) = 0.0;
    }

    is = 1;
    four1_FFTW(spec, ntime, is);

    //Output
    //ofstream cimout;
    //cimout.open ("test_CS(t).dat");

    //Prepare for passing to &vectors

    long int k = 0;  //will go from 0 to ntime-1
    //CRe.resize((int) (ntime));
    CIm.resize((int) (ntime));

    //w-Frame from [-ntime/2.0 -1 , ntime/2.0 -1]
    //for ntime = 8192*32 -> [-51471.5 , 51471.5]  which is maybe a bit too rough for current problems

    for(int i = ntime/2.0 -1; i>=1;i = i-1){
            //cimout<<-twopi*i/te<<" "<<spec(2*(int)ntime-2*i+1)*delta<<" "<<spec(2*(int)ntime-2*i+2)*delta<<" "<<endl;
            //CRe(k) = spec[2*(int)ntime-2*i+1]*delta;
            CIm(k) = spec(2*(int)ntime-2*i+2)*delta;
            k++;
    }

    for(int i = 0; i<=ntime/2-1;i++){
            //cimout<<twopi*i/te<<" "<<spec(2 * i +1) * delta<<" "<<spec(2 * i +2) *delta<<" "<<endl;
            //CRe(k) = spec[2 * i +1]*delta;
            CIm(k) = spec(2 * i +2)*delta;
            k++;
    }
    //cimout.close();
}


void PM_0_Four(int M, Vektor<double>& DMReals,Vektor<double> GReal, Vektor<double> GImag, double te, double ntime, Matrix<double> evecs, Vektor<double> evals, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat){

    double delta = te / (ntime -1);
    double om,spe,spe1;
    double twopi = 3.141592 * 2.0;
    double time = 0;
    double TauP = 550.0;

    cout << "TauP : " << TauP << endl;

    //Initialize Parameter
	Vektor<double> GReals;
	Vektor<double> GImags;

	GReals = GReal * y_MNKL(M, M, M, M, evecs, parameter_vec, dipolemat);
	GImags = GImag * y_MNKL(M, M, M, M, evecs, parameter_vec, dipolemat);

	double GM0 = GReals(0);
	double GM0i = GImags(0);

	cout << "Greals Size: " << GReals.size() << "   ntime - 1: " << ntime - 1 << endl;

    ublas::vector<double> spec(2 * (int)ntime + 1);

    //cout<<spec.size()<<endl;

    printf("Initializing Inputs\n");

    for(double i = 1; i <= ntime / 2.0 -1; i++){

//    	om = i * delta;

    	time = i * delta;

    	spe = exp(GReals(i) - GM0 - time * time / (4 * TauP * TauP)) * cos(GImags(i));
    	spe1 = exp(GReals(i) - GM0 - time * time / (4 * TauP * TauP)) * sin(GImags(i));

        spec(2 * i + 1) = spe;
        spec(2 * i + 2) = spe1;


        spe = exp(GReals(ntime - i) - GM0 - time * time / (4 * TauP * TauP)) * cos(GImags(ntime - i));
		spe1 = exp(GReals(ntime - i) - GM0 - time * time / (4 * TauP * TauP)) * sin(GImags(ntime - i));

		spec(2*(int)ntime - 2*i+1) = spe;
		spec(2*(int)ntime - 2*i+2) = spe1;



//    	spe = exp(-om*om/2.0);
//        spe1 = spe;
//
//        spec(2 * i + 1) = spe;
//        spec(2 * i + 2) = 0.0;
//
//        spec(2*(int)ntime - 2*i+1) = spe1;
//        spec(2*(int)ntime - 2*i+2) = 0.;
    }

    //spec(1) = 1.0 / (0.5*twopi);
    spec(1) = 1;
    spec(2) = 0;
    spec((int)ntime +1) = 0;
    spec((int)ntime +2) = 0;

    /*Nur für Testzwecke gebe ich hier die Funktionen aus, also die Originalfunktion und das was ich in den spec() Vektor gespeichert habe*/

//    ofstream test,test1;
//    test.open("Origfunktion.out");
//    test1.open("Sampledorigfunktion.out");
//
//    for(double i = 0; i < ntime ; i = i + 0.5){
//        test<<i<<" "<<exp(-i*i/2.0) <<endl;
//    }
//    for(double i = 0; i <ntime; i++){
//        test1<<twopi * i / delta / ntime<<" "<<spec(2*i+1)<<endl;
//    }
//
//    test1.close();
//    test.close();


    int is = -1;
    double corr_real, corr_imag;

    four1_FFTW(spec, ntime, is);

    printf("FFT Done! Doing write-outs!\n");

    ofstream FourTrans;
    FourTrans.open ("FourierTransformation.out");

    ublas::vector<double> corr(2 * (int)ntime + 1);

    for(double i = 1; i <= ntime; i++){

    	corr_real = spec(2*i-1) * delta;
    	corr_imag = spec(2*i) * delta;

//        corr_real = spec(2*i-1) / ntime * twopi / delta;
//        corr_imag = spec(2*i) / ntime * twopi / delta;
        //FourTrans <<setprecision(15)<< (i-1)*delta  <<" "<<corr_real<<" "<<corr_imag<<" "<<exp(-(i-1)*delta)<<endl;
//        FourTrans <<setprecision(15)<< (i-1)*delta  <<" "<<corr_real<<" "<<corr_imag<<" "<<sqrt(twopi)*exp(-((i-1)*delta)*((i-1)*delta)/2.0)<<endl;

        FourTrans <<setprecision(15)<< twopi*(i-1)/te  <<" "<<corr_real<<" "<<corr_imag<<" "<<sqrt(twopi)*exp(-twopi*(i-1)/te*twopi*(i-1)/te/2.0)<<endl;


        corr(2*i-1) = corr_real;
        corr(2*i) = corr_imag;
    }
    FourTrans.close();

    /*Die 2 Schleifen sind ein überbleibsel aus dem alten Code*/
    for(double i = 1; i<=ntime/2.0 ; i++){
            spec(2*i-1) = corr(2*i-1);
            spec(2*i) = corr(2*i);
    }

    for(double i = ntime/2.0 +1; i<=ntime; i++){
            spec(2*i-1) = corr(2*i-1);
            spec(2*i) = corr(2*i);
    }

//    printf("Doing Backwards FFT-Transformation!\n");
//
//    is = 1;
//    four1_FFTW(spec, ntime, is);
//
//    ofstream InvFourTrans;
//    InvFourTrans.open ("InvFourierTransformation.out");
//
//    for(double i = ntime/2.0 -1; i>=1;i = i-1){
//            //InvFourTrans<<-twopi*i/te<<" "<<spec(2*(int)ntime-2*i+1)*delta/twopi<<" "<<spec(2*(int)ntime-2*i+2)*delta<<" "<<1.0/(1.0 + twopi*i/te*twopi*i/te)*1.0 / (0.5*twopi)<<endl;
//            InvFourTrans<<-twopi*i/te<<" "<<spec(2*(int)ntime-2*i+1)*delta/twopi<<" "<<spec(2*(int)ntime-2*i+2)*delta/twopi<<" "<<exp(-twopi*i/te*twopi*i/te/2.0)<<endl;
//    }
//
//    for(double i = 0; i<=ntime/2-1;i++){
//
//            //InvFourTrans<<twopi*i/te<<" "<<spec(2 * i +1) * delta/twopi<<" "<<spec(2 * i +2) *delta<<" "<<1.0/(1.0 + twopi*i/te*twopi*i/te)*1.0 / (0.5*twopi)<<endl;
//            InvFourTrans<<twopi*i/te<<" "<<spec(2 * i +1) * delta/twopi<<" "<<spec(2 * i +2) *delta/twopi<<" "<<exp(-twopi*i/te*twopi*i/te/2.0)<<endl;
//    }
//
//    InvFourTrans.close();

    printf("Done!\n");
}


