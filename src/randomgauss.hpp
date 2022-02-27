#ifndef RANDOMGAUSS_HPP_INCLUDED
#define RANDOMGAUSS_HPP_INCLUDED

#include <bits/stdc++.h>
#include "ChlClass.hpp"
#include "random.hpp"
#include "Utilities.hpp"
#include "excmatrix.hpp"
#include "absorptionCD.hpp"
#include "correlation.hpp"
#include <omp.h>
#include "excdomains.hpp"
#include <stdio.h>
#include <fftw3.h>
#include "LU_Decomp.hpp"
#include "lineshape.hpp"

using namespace std;

double NormalRandom(double mean,double sigma){
 typedef boost::normal_distribution<double> NormalDistribution;
 typedef boost::mt19937 RandomGenerator;
 typedef boost::variate_generator<RandomGenerator&, \
                         NormalDistribution> GaussianGenerator;

  /* Initiate Random Number generator with seed 47382917263548729976411 */
  static RandomGenerator rng(static_cast<unsigned> (47382917263));

  /* Choose Normal Distribution */
  NormalDistribution gaussian_dist(mean, sigma);

  /* Create a Gaussian Random Number generator
   *  by binding with previously defined
   *  normal distribution object
   */
  GaussianGenerator generator(rng, gaussian_dist);

  // sample from the distribution
  return generator();
}

void RandomizeSiteEnergies(Vektor<double>& SiteEnergy, double FWHM){
    for(int i = 0; i < (int)SiteEnergy.size(); i++){
        SiteEnergy(i) = NormalRandom(SiteEnergy(i),FWHM);
    }
}

void RandomizeSiteEnergies_DiffBroadenings(Vektor<double>& SiteEnergy, Vektor<double>& Inh_Broads){

	//-----Static Disorder Site Energies with individual Inh. Broadenings for each Pigment-----
	for(int i = 0; i < (int)SiteEnergy.size(); i++){
        SiteEnergy(i) = NormalRandom(SiteEnergy(i), Inh_Broads(i));
    }
}

void Energytransfer_Redfield_GF(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, Vektor< ChlAtoms > Chl_Vec, Vektor<double>& Huang_Fac, Vektor<double>& Inh_Broads){
	//-----Energy Transfer using a Combination of Redfield and Generalized Förster Transfer with Exciton Domains-----

	//-----Initializations-----
	Vektor< Vektor<double> > Flu_sum ((long)parameter_vec(3));
	Vektor< Vektor< Vektor<double> >> Flu_sum_domains;
	double ntime = 8192.0 * 4.0;

	Vektor<double> GReals;
	Vektor<double> GImags;
	FT_Gt_WZ(GReals, GImags, 8192.0 / 2048.0, ntime, parameter_vec);

	Vektor<Matrix<double> > excmatrix_domains;
	Vektor<Vektor<double> > domains;

	//-----Ausgliederung von CIm & G(t)-----
	Vektor<double> CIm;
	correlation(parameter_vec, CIm);

	Vektor<Vektor<double> > coupling_vec;
	buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

	cout << Huang_Fac << endl;
	Vektor< Vektor<double> > Huang_Domains;
	Build_Huang_Domains(Huang_Fac, domains, Huang_Domains);

	Get_Coupling_Vec(coupling_vec);
	cout << domains << endl;

//	ofstream Jw;
//	Jw.open("J_w.out");
//	for(double w = 0; w < 100.0; w++){
//		Jw << w << " " << J_WZ(w) << endl;
//	}
//	Jw.close();
//
//
//	ofstream CReo;
//	CReo.open("CRe.out");
//	for(double w = 0; w < 100.0; w++){
//		CReo << w << " " << C_Re_WZ( w , parameter_vec(0) ) << endl;
//	}
//	CReo.close();
//	exit(1);
//	cout <<"Huang Domains: " << Huang_Domains << endl;

	printf("Initialization for Energy Transfer Done! Doing Static-Disorder Loop now.\n");

	long wrong_pop = 0;

	//-----Static-Disorder-----
	long ii = 0;
	#pragma omp parallel for private(ii) shared(Flu_sum, Flu_sum_domains)
	for(ii = 0; ii < (long)parameter_vec(3); ii++){

		if(ii % 50 == 0){
			cout << "Step #: " << ii << endl;
		}

		Vektor<double> sitevec_random;
		sitevec_random = sitevec;
//		vector<double> v = {12370.3,12460.6,12405,12381.4,12482.2,12592.9,12441.8,12792.4,12548.4,12394.5,12173.7,12344.7,12535.4,12476.2,12523.9,12674,12472.3,12474.9,12099.1,12378,12551.9,12495.6,12425.2,12633.7};
//		for(int i = 0; i < v.size(); i++){
//			sitevec_random(i) = v[i];
//		}

		//		cout << "sitevec: " << sitevec_random << endl;

		//-----The boost Normaldistribution has the inputform (MEAN-Value, Sigma)-----
		//-----parameter_vec(2)/2.35482 contains FWHM = 170cm^-1-----
		#pragma omp critical
		{
			RandomizeSiteEnergies(sitevec_random, parameter_vec(2) / 2.35482);
		}
		addsiteenergies(excmatrix_domains, domains, sitevec_random);

		Vektor< Matrix<double> > evecs_vec((long)excmatrix_domains.size());
		Vektor< Vektor<double> > evals_vec((long)excmatrix_domains.size());

//		cout << "Excmatrix_Domains: " << excmatrix_domains << endl;

		Vektor< Vektor< Vektor< Vektor<double> > > > dipolemat_domains ( (long)excmatrix_domains.size() );

		for(long d = 0; d < (long)excmatrix_domains.size(); d++){

			diagonalizeexcmatrix_domains(evecs_vec(d), evals_vec(d), chlname_vec, parameter_vec, excmatrix_domains(d));

			dipolemat_domains(d) = dipolemat;

			dipolemat_domains(d)(0).resize(domains(d).size());
			dipolemat_domains(d)(1).resize(domains(d).size());

			for(long u = 0; u < (long)dipolemat_domains(d)(0).size(); u++){

				dipolemat_domains(d)(0)(u) = dipolemat(0)(domains(d)(u)-1);
				dipolemat_domains(d)(1)(u) = dipolemat(1)(domains(d)(u)-1);
			}

		}

//		//-----Rescale the Pheos to |mu_Pheo| = 3.5 Debye-----
		dipolemat_domains(0)(0)(4) = dipolemat_domains(0)(0)(4) * 3.5 / 4.58;
		dipolemat_domains(0)(0)(5) = dipolemat_domains(0)(0)(5) * 3.5 / 4.58;

		//-----Initial Population PM(0)-----
		Vektor<double> P0;
		Vektor<double> DMRealsP0;

		for (long M = 0; M < (long) sitevec.size(); M++) {

			//-----Set initial Population to  1 / Npig, because sum_m |c_m^(M)|^2 = 1 -----
			push_back(P0, 1.0 / sitevec.size());
		}

		//-----Delay Time conversion fs -> cm-----
//		double delay = parameter_vec(7) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);

		//-----Build Kinetic Matrix A-----
		Matrix<double> A_kinetic((long)sitevec.size(), (long)sitevec.size());
		for (long M = 0; M < (long) A_kinetic.size1(); M++) {
			for (long N = 0; N < (long) A_kinetic.size2(); N++) {

				A_kinetic(M, N) = 0;
			}
		}

		long dom_count = 0;
		long dom_count_2 = 0;

		Vektor< Vektor <double> > A_Kinetic_Vektor;
		A_Kinetic_Vektor.resize((long)sitevec.size());
		for(long i = 0; i < (long)A_Kinetic_Vektor.size(); i++){
			A_Kinetic_Vektor(i).resize((long)sitevec.size());
		}

		//-----Build the Kinetic Matrix A with 2 Rates-----
		Vektor<double> DMReals;
		fftw_complex *in, *out;
		fftw_plan p;

		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

		for(long j = 0; j < evals_vec.size(); j++){
			reverse(evals_vec(j).begin(), evals_vec(j).end());
		}

//		Vektor<double> a(3);
//		a(0) = 0;
//		a(1) = 5;
//		a(2) = 10;
//
//		reverse(a.begin(), a.end());
//		cout << a << endl;

//		for(long j = 0; j < evecs_vec.size(); j++){
//
//			for(int M = 0; M < evecs_vec(j).size1(); M++){
//
//				Vektor<double>
//
//				for(int k = 0; k < evecs_vec(j).size2(); k++){
//
//				}
//			}
//		}
//
		for(long j = 0; j < evecs_vec.size(); j++){
			Matrix<double> mat_ordered;
			mat_ordered = evecs_vec(j);

			long s1 = evecs_vec(j).size1();
			long s2 = evecs_vec(j).size2();

			for(long k = 0; k < s1; k++){
				for(long l = 0; l < s2; l++){

//					mat_ordered(s1 - 1 - k, s2 - 1 -l) = evecs_vec(j)(k, l);
					mat_ordered(s1 - 1 - k, l) = evecs_vec(j)(k, l);
				}
			}

			evecs_vec(j) = mat_ordered;
//			cout << mat_ordered << endl;
		}

//		cout << evals_vec(0) << endl;
//		cout << evecs_vec(0) << endl;
//
//		cout << evecs_vec(0)(0,7) << endl;
//
//		cout << "done" << endl;
//		exit(1);

		//-----Vergleich bis hierher-----
		for(long a = 0; a < (long)domains.size(); a++){

			dom_count_2 = 0;

			for(long dd = 0; dd < (long)domains.size(); dd++){

				for(int M = 0; M < (int)domains(a).size(); M++){
					for(int N = 0; N <(int)domains(dd).size(); N++){
						//					cout << "a: " << a << " d: " << dd << " M: " << M << " N: " << N << " dom_count: " << dom_count << " dom_count_2: " << dom_count_2 << endl;
						if( a == dd){
							//-----Pigments are in the same domain-----

							if(M != N){
//								A_kinetic(dom_count + N, dom_count_2 + M) = -2.0 * y_MN_Huang(M, N, evecs_vec(a), Huang_Domains(a)) * C_Re_WZ( evals_vec(a)(M) - evals_vec(a)(N) , parameter_vec(0) );
//								A_kinetic(dom_count + N, dom_count_2 + M) = C_Re_WZ( evals_vec(a)(M) - evals_vec(a)(N) , parameter_vec(0) );
//								A_kinetic(dom_count + M, dom_count_2 + N) = C_Re_WZ( evals_vec(a)(M) - evals_vec(a)(N) , parameter_vec(0) );

								double rate = -2.0 * y_MN_Huang(M, N, evecs_vec(a), Huang_Domains(a)) * C_Re_WZ( evals_vec(a)(M) - evals_vec(a)(N) , parameter_vec(0) );
								if(abs(rate) < parameter_vec(8))
									rate = 0.0;

								A_kinetic(dom_count + M, dom_count_2 + N) = rate;

								//							A_Kinetic_Vektor(dom_count + N)(dom_count_2 + M)= -2.0 * y_MN_Huang(M, N, evecs_vec(a), Huang_Domains(a)) * C_Re_WZ( evals_vec(a)(M) - evals_vec(a)(N) , parameter_vec(0) );
//								cout << evecs_vec(a) << endl;
//								cout << evals_vec(a) << endl;
//								exit(1);
							}
						}
						if(a != dd){
							//-----Pigments are not in the same domain-----
							//-----Calculate the Foerster Rate-----

							double VMN = 0;
							double integ_value = 0;

//							integ_value = GF_kMK_WZ(N, M, dd, a, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
//							VMN = Foerster_Couplings(N, M, evecs_vec(dd), evecs_vec(a), coupling_vec, domains, dd, a);
//							A_kinetic(dom_count + M, dom_count_2 + N) = -2.0 * 3.141592 * VMN * VMN * integ_value;

							double wM0 = wM0_WZ_Huang(M, evals_vec(a), evecs_vec(a), parameter_vec, CIm, dipolemat, Huang_Domains(a));
							double wN0 = wM0_WZ_Huang(N, evals_vec(dd), evecs_vec(dd), parameter_vec, CIm, dipolemat, Huang_Domains(dd));

							if(wM0 - wN0 >= 0){
								integ_value = GF_kMK_WZ(M, N, a, dd, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
								VMN = Foerster_Couplings(M, N, evecs_vec(a), evecs_vec(dd), coupling_vec, domains, a, dd);
								double rate = -2.0 * 3.141592 * VMN * VMN * integ_value;

								//-----care about numeric stability-----
								if(abs(rate) < parameter_vec(8))
									rate = 0.0;

								A_kinetic(dom_count + M, dom_count_2 + N) = rate;
//								A_kinetic(dom_count + M, dom_count_2 + N) = integ_value;
							}
							else{
								double DB_fac = exp( (6.582 * pow(10.0,-16) * 1.88496 * pow(10.0,11) * (wM0-wN0)) / (8.617 * pow(10,-5) * parameter_vec(0) ) );
								integ_value = GF_kMK_WZ(N, M, dd, a, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
								VMN = Foerster_Couplings(N, M, evecs_vec(dd), evecs_vec(a), coupling_vec, domains, dd, a);
								double rate = -2.0 * 3.141592 * VMN * VMN * integ_value * DB_fac;

								//-----care about numeric stability-----
								if(abs(rate) < parameter_vec(8))
									rate = 0.0;

								A_kinetic(dom_count + M, dom_count_2 + N) = rate;
//								A_kinetic(dom_count + M, dom_count_2 + N) = integ_value * DB_fac;
							}


//							cout << M << " " << N << " " << evals_vec(a)(M) << " " << evals_vec(dd)(N) << " " << dom_count << " " << dom_count_2 << " " << wM0 << " " << wK0 << endl;


							//A_Kinetic_Vektor(dom_count + M)(dom_count_2 + N) = -2.0 * 3.141592 * VMN * VMN * integ_value;
//							cout << "M: " << M << "  N: " << N << " dom_a: " << a << " dom_dd: " << dd << "  [Foerster-Rate]" << "  integ_value: " << 2.0 * 3.141592 * integ_value << " VMN = " << VMN << " A(" << dom_count + M << "," << dom_count_2 + N << ")" << endl;
							//if(M == N){
							//	A_kinetic(dom_count + N, dom_count_2 + M) = -200.0;
							//}
						}
					}
				}
				dom_count_2 = dom_count_2 + domains(dd).size();
			}
			dom_count = dom_count + domains(a).size();
		}

//		int dom_t = 0;
//		cout << endl;
//		for(int M = 0; M < 8; M++){
//			double wM0 = wM0_WZ_Huang(M, evals_vec(dom_t), evecs_vec(dom_t), parameter_vec, CIm, dipolemat, Huang_Domains(dom_t));
//			cout << M << " " << evals_vec(dom_t)(M) << " " << wM0 << endl;
//		}

//		cout << fixed << setprecision(10) << A_kinetic << endl;
//		exit(1);
		//-----Replace Upwards Transfer by Detailed Balance Rates-----
//		//-----Test Detailed Balance on 2 Single-Pigment Domains-----
//		cout << "a -> b" << endl;
//		int dom_a = 1;
//		int dom_b = 2;
//		int M = 0;
//		int N = 0;
//
//		double integ_value = 0;
//		double VMN = 0;
//
//		cout << "here" << endl;
//		double wM0 = wM0_WZ_Huang(N, evals_vec(dom_b), evecs_vec(dom_b), parameter_vec, CIm, dipolemat, Huang_Domains(dom_b));
//		double wK0 = wM0_WZ_Huang(M, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, CIm, dipolemat, Huang_Domains(dom_a));
//
//		integ_value = GF_kMK_WZ(N, M, dom_b, dom_a, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
//		VMN = Foerster_Couplings(N, M, evecs_vec(dom_b), evecs_vec(dom_a), coupling_vec, domains, dom_b, dom_a);
//		cout << integ_value << endl;
//		cout << M << " " << N << " " << integ_value << " " << VMN << endl;
//		cout << wM0 << " " << wK0 << endl;
//
//		cout << "b -> a" << endl;
//		dom_a = 1;
//		dom_b = 2;
//		M = 0;
//		N = 0;
//
//		VMN = 0;
//		cout << "here" << endl;
//		wM0 = wM0_WZ_Huang(M, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, CIm, dipolemat, Huang_Domains(dom_a));
//		wK0 = wM0_WZ_Huang(N, evals_vec(dom_b), evecs_vec(dom_b), parameter_vec, CIm, dipolemat, Huang_Domains(dom_b));
//
//		integ_value = GF_kMK_WZ(M, N, dom_a, dom_b, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
//		VMN = Foerster_Couplings(M, N, evecs_vec(dom_b), evecs_vec(dom_a), coupling_vec, domains, dom_a, dom_b);
//		cout << integ_value << endl;
//		cout << M << " " << N << " " << integ_value << " " << VMN << endl;
//		cout << wM0 << " " << wK0 << endl;
//
//		//-----Alex's Constants-----
//		/** hc in cm*eV */
//		/** Boltzmann Konst. in [eV/k] */
//		double const hc2 = 1.239850e-4;
//		double const Boltz = 0.8617e-4;
//
//		double DB_rate = exp( (6.582 * pow(10.0,-16) * 1.88496 * pow(10.0,11) * (wM0-wK0)) / (8.617 * pow(10,-5) * parameter_vec(0) ) );
//		double DB_fac = exp( (wM0-wK0) * hc2 /(Boltz * parameter_vec(0)) );
//
//		cout << DB_fac << "  " << DB_rate << endl;
//		cout << integ_value << "  " << integ_value / DB_rate << endl;
//
//		cout << evals_vec(dom_a) << "  " << evals_vec(dom_b) << endl;
//		exit(1);
//
//		dom_count = 0;
//		for(long dom_a = 0; dom_a < (long)domains.size(); dom_a++){
//
//			dom_count_2 = 0;
//
//			for(long dom_b = 0; dom_b < (long)domains.size(); dom_b++){
//
//				for(long M = 0; M < (long)domains(dom_a).size(); M++){
//					for(long K = 0; K <(long)domains(dom_b).size(); K++){
//						//					cout << "a: " << a << " d: " << dd << " M: " << M << " N: " << N << " dom_count: " << dom_count << " dom_count_2: " << dom_count_2 << endl;
//						if( dom_a == dom_b){
//							//-----Pigments are in the same domain-----
//							continue;
//						}
//						if(dom_a != dom_b){
//							//-----Pigments are not in the same domain-----
//							//-----Calculate the Foerster Rate-----
//
//							double VMN = 0;
//							double integ_value = 0;
//
//							//-----Replace Coupling by DB Coupling-----
//							if(A_kinetic(dom_count + M, dom_count_2 + K) == 0){
//
//								double wM0 = wM0_WZ_Huang(M, evals_vec(dom_a), evecs_vec(dom_a), parameter_vec, CIm, dipolemat, Huang_Domains(dom_a));
//								double wK0 = wM0_WZ_Huang(K, evals_vec(dom_b), evecs_vec(dom_b), parameter_vec, CIm, dipolemat, Huang_Domains(dom_b));
//
//								double DB_fac = exp( (6.582 * pow(10.0,-16) * 1.88496 * pow(10.0,11) * (wM0-wK0)) / (8.617 * pow(10,-5) * parameter_vec(0) ) );
//
//								integ_value = GF_kMK_WZ(M, K, dom_a, dom_b, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in, out, p);
//								VMN = Foerster_Couplings(K, M, evecs_vec(dom_b), evecs_vec(dom_a), coupling_vec, domains, dom_b, dom_a);
//
//								A_kinetic(dom_count + M, dom_count_2 + K) = -2.0 * 3.141592 * VMN * VMN * integ_value / DB_fac;
//
////								cout << M << " " << K << " " << dom_count << " " << dom_count_2 << " " << integ_value << " " << VMN << endl;
////								exit(1);
//							}
//						}
//					}
//				}
//				dom_count_2 = dom_count_2 + domains(dom_b).size();
//			}
//			dom_count = dom_count + domains(dom_a).size();
//		}

		fftw_destroy_plan(p);
		fftw_free(in);
		fftw_free(out);

		//-----Sum Up the Diagonals-----
		auto A_kinetic_fin = trans(A_kinetic);
		for(long M = 0; M < (long)A_kinetic.size2(); M++){
			double k_sum = 0.0;

			for(long N = 0; N < (long)A_kinetic.size1(); N++){

				if(M != N){
					//				cout << "M: " << M << "  N: " << N << " A: " << A_kinetic(M, N) << endl;
					k_sum += A_kinetic(M, N);
				}
			}
			A_kinetic_fin(M, M) = -k_sum;
		}

		A_kinetic = -A_kinetic_fin;

//		for(long M = 0; M < (long)A_kinetic.size1(); M++){
//			double k_sum = 0.0;
//
//			for(long N = 0; N < (long)A_kinetic.size2(); N++){
//
//				if(M == N){
//					for(long K = 0; K < A_kinetic.size1(); K++){
//						k_sum -= A_kinetic(M, K);
//					}
//				}
//				else{
//					k_sum = A_kinetic(N, M);
//				}
//				A_kinetic_fin(M, N) = -k_sum;
//			}
//		}
//		for(long N = 0; N < (long)A_kinetic.size2(); N++){
//
//			double k_sum = 0.0;
//
//			for(long M = 0; M < (long)A_kinetic.size1(); M++){
//
//				if(M != N){
//					//				cout << "M: " << M << "  N: " << N << " A: " << A_kinetic(M, N) << endl;
//					k_sum += A_kinetic(M, N);
//
//					A_kinetic(N, N) = -k_sum;
//				}
//			}
//		}
//		A_kinetic = A_kinetic_fin;
		//-----Add the sink for CT in Domain 1 (Reaction Center)-----
		//-----The rate reads: k_M->RP1 = |c^(M)_ChlD1|^2 k_intr with k_intr = 1 / 6ps-----
		//-----c^(M)_ChlD1 is the AccD1 in Shibate et. al. and the 3rd Pigment in the list-----
		for(long M = 0; M < (long)domains(0).size(); M++){

			A_kinetic(M, M) = A_kinetic(M, M) - evecs_vec(0)(M, 2) * evecs_vec(0)(M, 2) * 1.0 / (parameter_vec(10) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926) );
			//		cout << "Rate for k_M->RP1 = " << evecs_vec(0)(M, 2) * evecs_vec(0)(M, 2) * 1.0 / (6000.0 * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926) ) << endl;
		}

		//-----Add the Radiative Terms in all domains-----
		//-----k_rad  = -|mu_M|^2 / |mu_0|^2 * 1.0/tau_Chl with tau_Chl = 4ns-----

		for (long d = 0; d < (long) domains.size(); d++) {

			long dom_size = 0;
			for(long cur_dom = 0; cur_dom < d; cur_dom++){
				dom_size += domains(cur_dom).size();
			}

			for (long M = 0; M < (long) domains(d).size(); M++) {

				Vektor<double> trans_dipole_moment(3);

				for(long j = 0; j < 3; j++)
					trans_dipole_moment(j) = 0.0;

				for (long j = 0; j < (long) domains(d).size(); j++) {
					trans_dipole_moment += dipolemat_domains(d)(0)(j) * evecs_vec(d)(M, j);
				}

				A_kinetic(M + dom_size, M + dom_size) = A_kinetic(M + dom_size, M + dom_size) - inner_prod(trans_dipole_moment, trans_dipole_moment) / inner_prod(dipolemat_domains(d)(0)(0), dipolemat_domains(d)(0)(0)) * 1.0 / (parameter_vec(9) * 1e3 * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926));
			}
		}

		Matrix<double> evecs_kinetic;
		Vektor<double> evals_kinetic;
		diagonalize_KineticMatrix(A_kinetic, evecs_kinetic, evals_kinetic);

		//-----Set Up the C_i's -> it's a Vektor of Vektors-----
		Vektor< Vektor<double> > ci_evecs;
		ci_evecs.resize((long)evecs_kinetic.size1());

		//-----Resize & Copy-----
		for(long i = 0; i < (long) ci_evecs.size(); i++){
			ci_evecs(i).resize((long)evecs_kinetic.size2());
		}

		for(long i = 0; i < (long) evecs_kinetic.size1(); i++){
			for(long j = 0; j < (long) evecs_kinetic.size2(); j++){
				ci_evecs(i)(j) = evecs_kinetic(i, j);
			}
		}
		//-----Normalize to 1 & Print-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){

			ci_evecs(i) = ci_evecs(i) / norm_2(ci_evecs(i));
			//    	cout << ci_evecs(i) / norm_2(ci_evecs(i)) << endl;
		}
		//-----Copy back the normalized Vektors to initial evecs_kinetic Matrix-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){
			for(long j = 0; j < (long)ci_evecs.size(); j++){

				evecs_kinetic(i, j) = ci_evecs(i)(j);
			}
		}
//		cout << "Normierte Kinetic Matrix: " << evecs_kinetic << endl;

		//-----Get Initial conditions d_i from P(0)-----
		//-----Solve the system of linear equations-----
		Matrix<double> c_matrix;
		c_matrix = trans(evecs_kinetic);

		Vektor<double> d;
		d = solve_lineq(c_matrix, P0);

//		cout << "Kinetic Matrix: " << A_kinetic << endl;
//		cout << "Initidal Conditions d_i: " << d << endl;
//		cout << "starting from Vektor P(0) = " << P0 << endl;
//		exit(1);
		//-----Output P(t)-----
		double sum = 0;
//		ofstream Pt;
//		Pt.open("PT_Vergleich_Problem.out");
//
//		double t_max = 1500;
//		double dt = 0.05;
//		int num_steps = (int)(t_max / dt);
//		Vektor< Vektor<double> > population(num_steps);
//
//		for(int j = 0; j < (int)population.size(); j++){
//			population(j).resize(evals_kinetic.size());
//			for(int k = 0; k < (int)population(j).size(); k++){
//				population(j)(k) = 0;
//			}
//		}
//
//		for(long k = 0; k < num_steps; k++){
//			Pt << k * dt << " ";
//			for(int l = 0; l < d.size(); l++){
//				sum = 0;
//				for(int j = 0; j < d.size(); j++){
//
//					sum = sum + ci_evecs(j)(l) * d(j) * 1.0 * exp(k * dt * evals_kinetic(j));
//				}
//				population(k)(l) += sum;
//				Pt << sum << " ";
//			}
//			Pt << endl;
//		}
//		Pt.close();
//
//		cout << "Population" << endl;
//
//		for(long k = 0; k < 100; k++){
//			int l = 16;
//			sum = 0;
//				for(int j = 0; j < d.size(); j++){
//					sum = sum + ci_evecs(j)(l) * d(j) * 1.0 * exp(k * dt * evals_kinetic(j));
//				}
//				cout << sum << endl;
//		}
//
//		cout << "A_kinetic: " << endl;
//		cout << A_kinetic << endl;
//
//		cout << evals_kinetic << endl;
//		cout << evecs_kinetic << endl;
//
//		cout << "Initial Condition" << endl;
//		cout << d << endl;
//
//		exit(1);
		//-----Integrate P(t) Vektor analytically -----
		Vektor<double> Pt_int((long)evals_kinetic.size());

		for (long i = 0; i < (long) Pt_int.size(); i++) {
			Pt_int(i) = 0;
		}

		for (long i = 0; i < (long) ci_evecs.size(); i++) {
			Pt_int = Pt_int + ci_evecs(i) * d(i) * 1.0 / (evals_kinetic(i) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926));
		}

		Pt_int = -Pt_int;

		for(int i = 0; i < Pt_int.size(); i++){
			if(Pt_int(i) < 0){
				wrong_pop++;
			}
		}

//		cout << Pt_int << endl;
//		cout << evals_kinetic << endl;
//		cout << fixed << setprecision(25) << A_kinetic << endl;
//		cout << setprecision(5) << sitevec_random << endl;
//		exit(1);
		//-----Fluorescence I(w) -----
		fftw_complex *in_, *out_;
		fftw_plan p_;

		in_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
		out_ = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

		Flu_sum_domains.resize(domains.size());
		for(long d = 0; d < (long) Flu_sum_domains.size(); d++){
			Flu_sum_domains(d).resize((long)parameter_vec(3));
		}

		for (long d = 0; d < (long) domains.size(); d++) {

			Vektor<double> DMReals;
			Vektor<double> Flu;

			long dom_size = 0;
			for(long cur_dom = 0; cur_dom < d; cur_dom++){
				dom_size += domains(cur_dom).size();
			}

			for (long M = 0; M < (long) domains(d).size(); M++) {

				Vektor<double> trans_dipole_moment(3);

				for(long j = 0; j < 3; j++)
					trans_dipole_moment(j) = 0.0;

				for (long j = 0; j < (long) domains(d).size(); j++) {
					trans_dipole_moment += dipolemat_domains(d)(0)(j) * evecs_vec(d)(M, j);
				}

				DM_Emission_WZ_Huang(M, d, DMReals, GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, CIm, evecs_vec, evals_vec, parameter_vec, dipolemat, Huang_Domains, in_, out_, p_);
				Flu.resize(DMReals.size());
				Flu = Flu + DMReals * inner_prod(trans_dipole_moment, trans_dipole_moment) * Pt_int(M + dom_size);

				if(Pt_int(M + dom_size) < 0){
					cout << "M: " << M << " d: " << d << " ii: " << ii << endl;
					cout << Pt_int << endl;
					cout << sitevec_random << endl;
				}

			}
			#pragma omp critical
			{
				Flu_sum(ii).resize(Flu.size());
				Flu_sum(ii) = Flu_sum(ii) + Flu;

				Flu_sum_domains(d)(ii).resize(Flu.size());
				Flu_sum_domains(d)(ii) = Flu_sum_domains(d)(ii) + Flu;
			}
		}

		fftw_destroy_plan(p_);
		fftw_free(in_);
		fftw_free(out_);

	}

	Vektor<double> Flu_sum_total;
	Flu_sum_total.resize(Flu_sum(0).size());

	for(long i = 0; i < (long)Flu_sum.size(); i++){
		Flu_sum_total = Flu_sum_total + Flu_sum(i);
	}

	Flu_sum_total = Flu_sum_total / parameter_vec(3);

	//-----Output to File-----
	cout << "Calculations done, writing out files!" << endl;

	FILE * Fludisorder;
	string FluFileName = "Flu_Dis.out";

	Fludisorder = fopen ( FluFileName.c_str() , "wb");
	for (long u = 0; u < (long) Flu_sum_total.size(); u++) {
		fprintf(Fludisorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), Flu_sum_total(u) );
	}
	fclose (Fludisorder);

	cout<<"Writing files done!"<<endl;

	if(wrong_pop > 0){
		cout << "Check Convergence of Balancing Algorithm for EV Calculation!!" << endl;
		cout << "There have been " << wrong_pop << " runs where this might be a problem!" << endl;
	}

	//-----Domain Contributions-----
	cout << "Writing out domain-contribution files!" << endl;
	cout << Flu_sum_domains.size() << endl;

	Vektor<Vektor<double>> Flu_sum_total_domains(Flu_sum_domains.size());
	for(long d = 0; d < Flu_sum_total_domains.size(); d++){
		Flu_sum_total_domains(d).resize(Flu_sum_domains(d)(0).size());
//		cout << "Flu_sum_domains(d)(0).size() " << Flu_sum_domains(d)(0).size() << endl;
	}

	for(long d = 0; d < Flu_sum_total_domains.size(); d++){
		for(long i = 0; i < (long)Flu_sum_domains(d).size(); i++){
			Flu_sum_total_domains(d) = Flu_sum_total_domains(d) + Flu_sum_domains(d)(i);
		}
	}
	for(long d = 0; d < Flu_sum_total_domains.size(); d++){
		Flu_sum_total_domains(d) = Flu_sum_total_domains(d) / parameter_vec(3);
	}
	//-----Output to File-----
//	cout << "Writing out domain-contribution files!" << endl;

	for(long d = 0; d < (long)Flu_sum_total_domains.size(); d++){
		FILE * Fludisorder_d;
		string FluFileName_d = "Flu_Dis_d" + to_string(d) + ".out";

		Fludisorder_d = fopen ( FluFileName_d.c_str() , "wb");
		for (long u = 0; u < (long) Flu_sum_total_domains(d).size(); u++) {
			fprintf(Fludisorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), Flu_sum_total_domains(d)(u) );
		}
		fclose (Fludisorder_d);
	}

	cout<<"Writing domain-contribution files done!"<<endl;


//	//-----Test Initial Conditions-----
//	Vektor<double> testP0(4);
//	for(int i = 0; i < 4; i++){
//		testP0(i) = 0;
//	}
//	for(int i = 0; i < 4; i++){
//		testP0(i) += d(0) * c_matrix(i, 0) + d(1) * c_matrix(i, 1) + d(2) * c_matrix(i, 2) + d(3) * c_matrix(i, 3);
//	}
//	cout << "Testing Initial Conditions: " << testP0 << endl;
//	Matrix<double> A(4,4);
//	A(0,0) = 1; A(0,1) = 1; A(0,2) = 1; A(0,3) = 10;
//	A(1,0) = 1; A(1,1) = 2; A(1,2) = 3; A(1,3) = 7;
//	A(2,0) = 1; A(2,1) = 5; A(2,2) = 6; A(2,3) = 4;
//	A(3,0) = 4; A(3,1) = 5; A(3,2) = 7; A(3,3) = 8;
//
//	Vektor<double> b(4);
//	b(0) = 1; b(1) = 2; b(2) = 3; b(3) = 4;
//
//	Vektor<double> x;
//	x = solve_lineq(A, b);
//	cout << "Solved x: " << x << endl;

}

double StaticDisorder_domains(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, double t_ang){
    //Static Disorder


    Vektor< Vektor<double> > alpha_sum ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
    Vektor< Vektor<double> > RotStrength_sum ((int)parameter_vec(3));


    Vektor< Matrix<double> > excmatrix_domains;
    Vektor< Vektor<double> > domains;

    //Ausgliederung von CIm & G(t)
    Vektor<double> CIm;
    correlation(parameter_vec, CIm);

    Vektor<double> GReals;
    Vektor<double> GImags;
    FT_Gt_WZ(GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,parameter_vec);    //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0

    Vektor< Vektor<double> > coupling_vec;
//    build_predefined_excitondomains(chlname_vec, excmatrix_domains, domains, coupling_vec);
    buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    //cout<<"Domainentry :  "<<domains(0)(0)<<endl;
    //cout<<endl<<endl<<"Dipolemat :  "<<dipolemat<<endl;
    //cout<<"Excmatrix_domains size :  "<<excmatrix_domains.size()<<endl;

    double ntime = 8192.0 * 4.0;

    printf("Initialization Done! Doing Static-Disorder Loop now.\n");

    int i = 0;
//	#pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
	#pragma omp parallel for
    for(i = 0; i < (int)parameter_vec(3); i++){ //parameter_vec(3) holds #nrand

        Vektor<double> CD;
        Vektor<double> alpha;
        Matrix<double> evecs;            //matrix&vector for Eigenvectors & Eigenvalues
        Vektor<double> evals;

        Vektor<double> RotStrength_vec;

        Vektor<double> sitevec_random;
        sitevec_random = sitevec;

        fftw_complex *in, *out;
        fftw_plan p;

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

        RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
//        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
        addsiteenergies(excmatrix_domains,domains,sitevec_random);

        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;

        for(int j = 0; j < (int)excmatrix_domains.size(); j++){  //(int)excmatrix_domains.size()

            Vektor< Vektor< Vektor<double> > > dipolemat_domains;
            dipolemat_domains = dipolemat;

            //cout<<"Building dipolemat_domains :  "<<endl;

            dipolemat_domains(0).resize(domains(j).size());
            dipolemat_domains(1).resize(domains(j).size());

            for(int u = 0; u < (int)dipolemat_domains(0).size(); u++){

                dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u)-1);
                dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u)-1);
            }


            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));

//            get_RotStrength(evecs, evals, dipolemat_domains, parameter_vec, RotStrength_vec, 14814.8);

            CDabsorption_domains(evecs, evals, dipolemat_domains, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);


			#pragma omp critical
            {
                CD_sum(i) = CD;
                alpha_sum(i) = alpha;
                RotStrength_sum(i) = RotStrength_vec;
            }
            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
        }

        sitevec_random = sitevec;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }

//    FILE * RotStrength_dis_out;
//    RotStrength_dis_out = fopen ("RotStrength.out", "a+");
//
//    fprintf(RotStrength_dis_out,"<R->dis = %lf at Angle %lf \n", dis_RotStrength(RotStrength_sum, -1), t_ang);
//    cout << "dis_RotStrength: " << dis_RotStrength(RotStrength_sum, -1) << endl;
//
//    fclose(RotStrength_dis_out);

    Vektor<double> CD_sum_total;
    Vektor<double> alpha_sum_total;

    CD_sum_total.resize(CD_sum(0).size());
    alpha_sum_total.resize(alpha_sum(0).size());

    for(int i = 0; i < (int)CD_sum.size(); i++){
        CD_sum_total = CD_sum_total + CD_sum(i);
        alpha_sum_total = alpha_sum_total + alpha_sum(i);
    }

    alpha_sum_total = alpha_sum_total / parameter_vec(3);
    CD_sum_total = CD_sum_total / parameter_vec(3);

    cout<<"Calculations done, writing out files!"<<endl;

    FILE * ODdisorder;
    FILE * CDdisorder;

    string ODFileName = "OD_Dis.out";
    string CDFileName = "CD_Dis.out";

    ODdisorder = fopen ( ODFileName.c_str() , "wb");
    CDdisorder = fopen ( CDFileName.c_str() , "wb");

    for(int u = 0; u < (int)CD_sum_total.size(); u++){

        fprintf(ODdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),alpha_sum_total(u));  //16
        fprintf(CDdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),CD_sum_total(u));     //16
    }

    fclose (ODdisorder);
    fclose (CDdisorder);

    cout<<"Writing files done!"<<endl;

    return 0;
}

double StaticDisorder_domains_Get_Contributions(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, double t_ang){
    //Static Disorder

    Vektor< Vektor<double> > alpha_sum_all ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_all ((int)parameter_vec(3));

    Vektor< Vektor<double> > alpha_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > alpha_sum_HE ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_HE ((int)parameter_vec(3));

    Vektor< Vektor<double> > RotStrength_sum ((int)parameter_vec(3));


    Vektor< Matrix<double> > excmatrix_domains;
    Vektor< Vektor<double> > domains;

    //Ausgliederung von CIm & G(t)
    Vektor<double> CIm;
    correlation(parameter_vec, CIm);

    Vektor<double> GReals;
    Vektor<double> GImags;
    FT_Gt_WZ(GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,parameter_vec);    //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0

    Vektor< Vektor<double> > coupling_vec;
//    build_predefined_excitondomains(chlname_vec, excmatrix_domains, domains, coupling_vec);
    buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    //cout<<"Domainentry :  "<<domains(0)(0)<<endl;
    //cout<<endl<<endl<<"Dipolemat :  "<<dipolemat<<endl;
    //cout<<"Excmatrix_domains size :  "<<excmatrix_domains.size()<<endl;

    double ntime = 8192.0 * 4.0;

    printf("Initialization Done! Doing Static-Disorder Loop now.\n");

    int i = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
    for(i = 0; i < parameter_vec(3); i++){ //parameter_vec(3) holds #nrand

        Vektor< Vektor<double> > CD(3);
        Vektor< Vektor<double> > alpha(3);
        Matrix<double> evecs;            //matrix&vector for Eigenvectors & Eigenvalues
        Vektor<double> evals;

        Vektor<double> RotStrength_vec;

        Vektor<double> sitevec_random;
        sitevec_random = sitevec;

        fftw_complex *in, *out;
        fftw_plan p;

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);
//        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
        addsiteenergies(excmatrix_domains,domains,sitevec_random);

        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;

        for(int j = 0; j < (int)excmatrix_domains.size(); j++){  //(int)excmatrix_domains.size()

            Vektor< Vektor< Vektor<double> > > dipolemat_domains;
            dipolemat_domains = dipolemat;

            //cout<<"Building dipolemat_domains :  "<<endl;

            dipolemat_domains(0).resize(domains(j).size());
            dipolemat_domains(1).resize(domains(j).size());

            for(int u = 0; u < (int)dipolemat_domains(0).size(); u++){

                dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u)-1);
                dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u)-1);
            }


            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));

//            get_RotStrength(evecs, evals, dipolemat_domains, parameter_vec, RotStrength_vec, 14814.8);

//            CDabsorption_domains(evecs, evals, dipolemat_domains, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);
            CDabsorption_domains_Get_Contributions(evecs, evals, dipolemat, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);


//            #pragma omp critical
            {
                CD_sum_all(i) = CD(0);
                alpha_sum_all(i) = alpha(0);

                CD_sum_LE(i) = CD(1);
                alpha_sum_LE(i) = alpha(1);

                CD_sum_HE(i) = CD(2);
                alpha_sum_HE(i) = alpha(2);

                RotStrength_sum(i) = RotStrength_vec;
            }
            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
        }

        sitevec_random = sitevec;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }

//    FILE * RotStrength_dis_out;
//    RotStrength_dis_out = fopen ("RotStrength.out", "a+");
//
//    fprintf(RotStrength_dis_out,"<R->dis = %lf at Angle %lf \n", dis_RotStrength(RotStrength_sum, -1), t_ang);
//    cout << "dis_RotStrength: " << dis_RotStrength(RotStrength_sum, -1) << endl;
//
//    fclose(RotStrength_dis_out);

    cout << "Loops done! Summing up before writing out..." << endl;

    Vektor<double> CD_sum_total;
    Vektor<double> CD_sum_total_LE;
    Vektor<double> CD_sum_total_HE;
    Vektor<double> alpha_sum_total;
    Vektor<double> alpha_sum_total_LE;
    Vektor<double> alpha_sum_total_HE;

    CD_sum_total.resize(CD_sum_all(0).size());
    alpha_sum_total.resize(alpha_sum_all(0).size());

    CD_sum_total_LE.resize(CD_sum_LE(0).size());
    alpha_sum_total_LE.resize(alpha_sum_LE(0).size());

    CD_sum_total_HE.resize(CD_sum_HE(0).size());
    alpha_sum_total_HE.resize(alpha_sum_HE(0).size());

    for(int i = 0; i < (int)CD_sum_all.size(); i++){

        CD_sum_total = CD_sum_total + CD_sum_all(i);
        alpha_sum_total = alpha_sum_total + alpha_sum_all(i);

        CD_sum_total_LE = CD_sum_total_LE + CD_sum_LE(i);
        alpha_sum_total_LE = alpha_sum_total_LE + alpha_sum_LE(i);

        CD_sum_total_HE = CD_sum_total_HE + CD_sum_HE(i);
        alpha_sum_total_HE = alpha_sum_total_HE + alpha_sum_HE(i);

    }


    alpha_sum_total = alpha_sum_total / parameter_vec(3);
    CD_sum_total = CD_sum_total / parameter_vec(3);

    alpha_sum_total_LE = alpha_sum_total_LE / parameter_vec(3);
    CD_sum_total_LE = CD_sum_total_LE / parameter_vec(3);

    alpha_sum_total_HE = alpha_sum_total_HE / parameter_vec(3);
    CD_sum_total_HE = CD_sum_total_HE / parameter_vec(3);

    cout<<"Calculations done, writing out files!"<<endl;

    FILE * ODdisorder;
    FILE * CDdisorder;

    string ODFileName = "OD_Dis.out";
    string CDFileName = "CD_Dis.out";

    ODdisorder = fopen ( ODFileName.c_str() , "wb");
    CDdisorder = fopen ( CDFileName.c_str() , "wb");

    for(int u = 0; u < (int)CD_sum_total.size(); u++){

		fprintf(ODdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), alpha_sum_total(u), alpha_sum_total_LE(u), alpha_sum_total_HE(u));  //16
		fprintf(CDdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), CD_sum_total(u), CD_sum_total_LE(u), CD_sum_total_HE(u));     //16
    }

    fclose (ODdisorder);
    fclose (CDdisorder);

    cout<<"Writing files done!"<<endl;

    return 0;
}

void Lin_Pump_Test(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat){
	//-----Linear Pump-Test with Static Disorder-----

	//-----Initializations-----

	Vektor<Vektor<double> > GB_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > SE_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > ESA_sum((int) parameter_vec(3));

	double ntime = 8192.0 * 4.0;

	Vektor<double> GReals;
	Vektor<double> GImags;
	FT_Gt_WZ(GReals, GImags, 8192.0 / 2048.0, ntime, parameter_vec);

	Vektor<Matrix<double> > excmatrix_domains;
	Vektor<Vektor<double> > domains;

	Vektor<double> CIm;
	correlation(parameter_vec, CIm);

	Vektor<Vektor<double> > coupling_vec;
	buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    printf("Initialization for Pump-Test Done! Doing Static-Disorder Loop now.\n");

    //-----Insert Static-Disorder-----
    int i = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
//    #pragma omp parallel for
    for(i = 0; i < parameter_vec(3); i++){

    	if(i % 50 == 0){
    		cout << "DisorderLoop at i = " << i << endl;
    	}

    	//-----Building 1-Exciton Matrices-----
    	Matrix<double> evecs;
    	Vektor<double> evals;

    	Vektor<double> sitevec_random;
    	sitevec_random = sitevec;

    	//-----The boost Normaldistribution has the inputform (MEAN-Value, Sigma)-----
    	//-----parameter_vec(2)/2.35482 contains FWHM = 170cm^-1-----

    	RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
    	addsiteenergies(excmatrix_domains, domains, sitevec_random);

    	diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(0));

    	//-----Building 2-Exciton Matrices-----
    	Matrix<double> Exc2Matrix;
    	Matrix<double> evecs2N;
    	Vektor<double> evals2N;
    	Matrix<double> Exc2_Map;

    	Build_2Exiton_Matrix(excmatrix_domains(0), Exc2Matrix);

    	diagonalizeexcmatrix_domains(evecs2N, evals2N, chlname_vec, parameter_vec, Exc2Matrix);

    	Build_2Exc_Evec_Map(excmatrix_domains(0), Exc2_Map);

    	//    cout <<"Exc2 Matrix : " << Exc2Matrix << endl;
    	//    cout << "Evals2Exc: " << evals2N << "   Evecs(0,0): " << evecs2N(0,0) << endl;
    	//    cout << "Evals2Exc size(): " << evals2N.size() << " evecsizes: " << evecs2N.size1() << "  " << evecs2N.size2() << endl;
    	//    //double y_2N2K (int N2, int K2, Matrix<double> evec2Exc, Matrix<double> Exc2_Map, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat){
    	//
    	//    cout << "y2N2N: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "w2N0 : " << w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat) << endl;
    	//
    	//    double w2N0 = w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat);
    	//
    	//    cout << "1-Exciton Evecs: " << evecs << endl;
    	//
    	//    cout << "y2N2NMM: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "y2N2NMM: " << y_2N2KMM(0, 1, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//
    	//    cout << "tau2N: " << tau2N_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//
    	//    cout << "G(t) Prefactor: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "G(t) Prefactor: " << y_MNKL(0, 0, 0, 0, evecs, parameter_vec, dipolemat) << endl;
    	//    cout << "G(t) Prefactor: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;

    	//-----Initialize the Transition Dipole Moments-----
    	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
    	ublas::zero_vector<double> zero(3);

    	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {

    		trans_dipole_moment(M) = zero;
    	}

    	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {
    		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {

    			trans_dipole_moment(M) += dipolemat(0)(j) * evecs(M, j);
    		}
    	}

    	Matrix< Vektor<double> > trans_dipole_moment_2N( (int)evals2N.size(), (int)evals.size() );

    	for (int M = 0; M < (int) evals.size(); M++) {

    		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

    			trans_dipole_moment_2N(N2, M) = zero;
    		}
    	}

    	for (int M = 0; M < (int) evals.size(); M++) {
    		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

    			for (int k = 0; k < (int) evals.size(); k++) {
    				for (int l = 0; l < (int) evals.size(); l++) {

    					if (l > k) {

    						trans_dipole_moment_2N(N2, M) += evecs2N(N2, Exc2_Map(k, l)) * ( dipolemat(0)(k) * evecs(M,l) + dipolemat(0)(l) * evecs(M,k) );
    						//						cout << "trans_dipole_moment_2N" << trans_dipole_moment_2N(N2, M) << " evecs2N " << evecs2N(N2, Exc2_Map(k, l)) << " mu_k: " <<  dipolemat(0)(k) * evecs(M,k) << " mu_l: " << dipolemat(0)(l) * evecs(M,l) << endl;
    					}

    				}
    			}
    		}

    	}

    	//-----Initial Population PM(0)-----
    	Vektor<double> P0;

    	//-----Simple Dimer-----
    	push_back(P0, 1.0);
    	push_back(P0, 0.0);
//    	cout << "evals & p0: " << endl;
//    	cout << evals << endl;
//    	cout << P0 << endl;

//    	double U_w = 0.0;
//    	Vektor<double> DMRealsP0;
//    	fftw_complex *in3, *out3;
//
//    	in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//    	out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//
//    	for (int M = 0; M < (int) evals.size(); M++) {
//
//    		U_w = PM0_WZ(M, DMRealsP0, GReals, GImags, 8192.0 / 2048.0, ntime, CIm, evecs, evals, parameter_vec, dipolemat, in3, out3, parameter_vec(6), parameter_vec(5));
//    		//push_back(P0, U_w * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) );
//    		push_back(P0, U_w * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) );
//    	}
//
//    	fftw_free(in3);
//    	fftw_free(out3);
//
//    	//-----DEBUG-----
//    	cout <<"DEBUG -- P0: " << P0 << " " << P0 / norm_2(P0) << endl;
//    	cout << "mu_0: " << inner_prod(trans_dipole_moment(0), trans_dipole_moment(0)) << endl;
//    	cout << "mu_1: " << inner_prod(trans_dipole_moment(1), trans_dipole_moment(1)) << endl;
//    	exit(1);
    	//
    	//	//-----DEBUG-----
    	//	cout <<"DEBUG -- Evals / Evecs: " << evals << endl << evecs << endl;
    	//
    	//	//-----DEBUG-----
    	//	cout <<"DEBUG -- Evals2N / Evecs2N: " << evals2N << endl << evecs2N << endl;
    	//
    	//	//-----DEBUG-----
    	//	cout <<"DEBUG -- Transition Dipoles: " << trans_dipole_moment << endl;
    	//
    	//	//-----DEBUG-----
    	//	cout <<"DEBUG -- Transition Dipoles 2N: " << trans_dipole_moment_2N << endl;

    	//-----Solve Rate Equation for Kinetic Matrix A-----

    	//-----Build Kinetic Matrix A-----
    	Matrix<double> A_kinetic((long)sitevec.size(), (long)sitevec.size());
    	for (long M = 0; M < (long) A_kinetic.size1(); M++) {
    		for (long N = 0; N < (long) A_kinetic.size2(); N++) {

    			A_kinetic(M, N) = 0;
    		}
    	}

    	cout << "evals(M) - evals(N): " << evals(0) - evals(1) << endl;
    	cout << -2.0 * y_MN(0, 1, evecs) * C_Re_WZ(evals(0) - evals(1), parameter_vec(0) ) << endl;


		for(int M = 0; M < (int)evals.size(); M++){
			for(int N = 0; N < (int)evals.size(); N++){

				if(M != N){
					double rate = -2.0 * y_MN(M, N, evecs) * C_Re_WZ(evals(M) - evals(N), parameter_vec(0) );
					A_kinetic(M, N) = rate;
				}
			}
		}

		auto A_kinetic_fin = trans(A_kinetic);

		cout << A_kinetic_fin << endl;

		//-----Sum Up the Diagonals-----

//		auto A_kinetic_fin = A_kinetic;
		for(long N = 0; N < (long)A_kinetic_fin.size1(); N++){
			double k_sum = 0.0;
			for(long M = 0; M < (long)A_kinetic_fin.size2(); M++){

					if(M != N){
						k_sum += A_kinetic_fin(M, N);
					}
				}
				A_kinetic_fin(N, N) = -k_sum;
		}

		A_kinetic = -A_kinetic_fin;

		cout << "A_kinetic : " << A_kinetic << endl;

		Matrix<double> evecs_kinetic;
		Vektor<double> evals_kinetic;

//		for(long i = 0; i < (long) evecs_kinetic.size1(); i++){
//			cout << "{";
//			for(long j = 0; j < (long) evecs_kinetic.size2(); j++){
//				cout << evecs_kinetic(i, j) << ",";
//			}
//			cout << "}" << endl;
//		}
		diagonalize_KineticMatrix(A_kinetic, evecs_kinetic, evals_kinetic);

		//-----Set Up the C_i's -> it's a Vektor of Vektors-----
		Vektor< Vektor<double> > ci_evecs;
		ci_evecs.resize((long)evecs_kinetic.size1());

		//-----Resize & Copy-----
		for(long i = 0; i < (long) ci_evecs.size(); i++){
			ci_evecs(i).resize((long)evecs_kinetic.size2());
		}

		for(long i = 0; i < (long) evecs_kinetic.size1(); i++){
			for(long j = 0; j < (long) evecs_kinetic.size2(); j++){
				ci_evecs(i)(j) = evecs_kinetic(i, j);
			}
		}
		//-----Normalize to 1 & Print-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){

			ci_evecs(i) = ci_evecs(i) / norm_2(ci_evecs(i));
			cout << ci_evecs(i) / norm_2(ci_evecs(i)) << endl;
		}
		//-----Copy back the normalized Vektors to initial evecs_kinetic Matrix-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){
			for(long j = 0; j < (long)ci_evecs.size(); j++){

				evecs_kinetic(i, j) = ci_evecs(i)(j);
			}
		}
		//		cout << "Normierte Kinetic Matrix: " << evecs_kinetic << endl;

		//-----Get Initial conditions d_i from P(0)-----
		//-----Solve the system of linear equations-----
		Matrix<double> c_matrix;
		c_matrix = trans(evecs_kinetic);

		Vektor<double> d;
		d = solve_lineq(c_matrix, P0);

		cout << "c_matrix: " << c_matrix << endl;
		cout << "Initidal Conditions d_i: " << d << endl;
		cout << "starting from Vektor P(0) = " << P0 << endl;
		cout << "evals_kinetic: " << evals_kinetic << endl;
		cout << "evecs_kinetic: " << evecs_kinetic << endl;

		//-----Output P(t)-----
//		double sum = 0;
//
//		ofstream Pt;
//		Pt.open("PT_Vergleich_Problem.out");
//
//		double t_max = 150;
//		double dt = 0.005;
//		int num_steps = (int)(t_max / dt);
//		Vektor< Vektor<double> > population(num_steps);
//
//		for(int k = 0; k < (int)population.size(); k++){
//			population(k).resize(evals.size());
//		}
//
//		for(long k = 0; k < num_steps; k++){
//			Pt << k * dt << " ";
//			for(int i = 0; i < (int)d.size(); i++){
//
//				population(k) += d(i) * ci_evecs(i) * exp(dt * k * evals_kinetic(i));
//			}
//			Pt << population(k) << endl;
//		}
//
//		Pt.close();

		//-----Delay Time conversion fs -> cm-----
		double delay = parameter_vec(7) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);
		Vektor<double> PM_delay(evals_kinetic.size());

		for(int i = 0; i < (int)d.size(); i++){
			PM_delay += d(i) * ci_evecs(i) * exp(delay * evals_kinetic(i));
		}

		cout << delay << endl;
		cout << PM_delay << endl;

//		//-----Equilibrium Populations-----
//
//		double t_inf = 50000 * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);
//		Vektor<double> PM_equ(evals_kinetic.size());
//		Vektor<double> PM_equ_2(evals_kinetic.size());
//
//		for(int i = 0; i < (int)d.size(); i++){
//			PM_equ += d(i) * ci_evecs(i) * exp(t_inf * evals_kinetic(i));
//		}
//
//		for(int i = 0; i < (int)d.size(); i++){
//			PM_equ_2 += d(i) * ci_evecs(i) * exp(t_inf * 2 * evals_kinetic(i));
//		}
//
//		//-----Check if equilibrium is reached-----
//		double error = 0;
//		for(int M = 0; M < (int)PM_equ.size(); M++){
//			error += abs(PM_equ(M) - PM_equ_2(M));
//		}
//
//		if(error > 1e-3){
//			cout << "DEBUG: Equlibration-Error: " << error << endl;
//			cout << "Check Equilibration-Time!! Exiting for now!" << endl;
//			exit(1);
//		}

		//-----Groundstate Bleaching-----
    	Vektor<double> GB;
    	fftw_complex *inGB, *outGB;

    	inGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
    	outGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

    	Lin_GSBleach(evecs, evals, dipolemat, parameter_vec, GB, CIm, GReals, GImags, inGB, outGB, P0, trans_dipole_moment);

    	//-----GB Output Stuff-----
    	FILE * ODdisorder;
    	string ODFileName = "GB_Test.out";
    	ODdisorder = fopen(ODFileName.c_str(), "wb");

    	for (int w = 0; w < 15000; w++) {
    		fprintf(ODdisorder, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , GB(w));
    	}
    	fclose(ODdisorder);

//    	for (int u = 0; u < (int) GB.size(); u++) {
//
//    		fprintf(ODdisorder, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))),GB(u));
//    	}
//    	fclose(ODdisorder);

    	//-----Stimulated Emission-----
    	Vektor<double> SE;
    	fftw_complex *inSE, *outSE;

    	inSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
    	outSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

    	Lin_StimEmission(evecs, evals, dipolemat, parameter_vec, SE, CIm, GReals, GImags, inSE, outSE, P0, trans_dipole_moment, PM_delay);

    	//-----SE Output Stuff-----
    	FILE * SEDis;
    	string SEFileName = "SE_Test.out";
    	SEDis = fopen(SEFileName.c_str(), "wb");

    	for (int w = 0; w < 15000; w++) {
    		fprintf(SEDis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , SE(w));
    	}

//    	for (int u = 0; u < (int) SE.size(); u++) {
//
//    		fprintf(SEDis, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), SE(u));
//    	}
    	fclose(SEDis);

//    	exit(1);
    	//-----Excited State Absorption-----
    	Vektor<double> ESA;
    	fftw_complex *inESA, *outESA;

    	inESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
    	outESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

    	Lin_ESA(evecs, evals, evecs2N, evals2N, Exc2_Map, dipolemat, parameter_vec, ESA, CIm, GReals, GImags, inESA, outESA, P0, trans_dipole_moment, trans_dipole_moment_2N, PM_delay);

    	//-----ESA Output Stuff-----
    	FILE * ESADis;
    	string ESAFileName = "ESA_Test.out";
    	ESADis = fopen(ESAFileName.c_str(), "wb");

    	for (int w = 0; w < 15000; w++) {
    		fprintf(ESADis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , ESA(w));
    	}

//    	for (int u = 0; u < (int) ESA.size(); u++) {
//
//    		fprintf(ESADis, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ESA(u));
//    	}
    	fclose(ESADis);

    	#pragma omp critical
    	{
    		GB_sum(i) = GB;
    		SE_sum(i) = SE;
    		ESA_sum(i) = ESA;
    	}

    	sitevec_random = sitevec;

    	fftw_free(inGB);
    	fftw_free(outGB);

    	fftw_free(inSE);
    	fftw_free(outSE);

    	fftw_free(inESA);
    	fftw_free(outESA);

    }

    Vektor<double> GB_sum_total;
	Vektor<double> SE_sum_total;
	Vektor<double> ESA_sum_total;
	Vektor<double> Pump_Test_total;

	GB_sum_total.resize(GB_sum(0).size());
	SE_sum_total.resize(SE_sum(0).size());
	ESA_sum_total.resize(ESA_sum(0).size());
	Pump_Test_total.resize(ESA_sum(0).size());

	for (int j = 0; j < (int) GB_sum.size(); j++) {
		GB_sum_total = GB_sum_total + GB_sum(j);
		SE_sum_total = SE_sum_total + SE_sum(j);
		ESA_sum_total = ESA_sum_total + ESA_sum(j);
	}

	GB_sum_total = GB_sum_total / parameter_vec(3);
	SE_sum_total = SE_sum_total / parameter_vec(3);
	ESA_sum_total = ESA_sum_total / parameter_vec(3);

	Pump_Test_total = GB_sum_total + SE_sum_total + ESA_sum_total;

	cout << "Calculations done, writing out files!" << endl;

	FILE * Pump_Test_Disorder;
	string Pump_Test_FileName = "Pump_Test_Dis_" + to_string((int)parameter_vec(7)) + ".out";

	Pump_Test_Disorder = fopen(Pump_Test_FileName.c_str(), "wb");

//	for (int u = 0; u < (int) GB_sum_total.size(); u++) {
//
//		fprintf(Pump_Test_Disorder, " %lf  %lf  %lf  %lf  %lf\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), Pump_Test_total(u), GB_sum_total(u), SE_sum_total(u), ESA_sum_total(u));  //16
//		//-----w, Pump-Test, GB, SE, ESA-----
//	}

	for (int w = 0; w < 15000; w++) {
		fprintf(Pump_Test_Disorder, " %lf  %lf  %lf  %lf  %lf\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , Pump_Test_total(w), GB_sum_total(w), SE_sum_total(w), ESA_sum_total(w));
	}

	fclose(Pump_Test_Disorder);
	cout << "Writing files done!" << endl;
}

void CD_Pump_Test(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, Matrix<double>& Exc2Shifts){
	//-----TRCD (Circular Pump-Test), Linear Pump-Test and CD / OD with Static Disorder-----

	//-----Initializations-----

	Vektor<Vektor<double> > GB_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > SE_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > ESA_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > CD_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > OD_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > CDGB_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > CDSE_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > CDESA_sum((int) parameter_vec(3));


	double ntime = 8192.0 * 4.0;

	Vektor<double> GReals;
	Vektor<double> GImags;
	FT_Gt_WZ(GReals, GImags, 8192.0 / 2048.0, ntime, parameter_vec);

	Vektor<Matrix<double> > excmatrix_domains;
	Vektor<Vektor<double> > domains;

	Vektor<double> CIm;
	correlation(parameter_vec, CIm);

	Vektor<Vektor<double> > coupling_vec;
	buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    printf("Initialization for Pump-Test Done! Doing Static-Disorder Loop now.\n");

    //-----Insert Static-Disorder-----
    int ii = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)

    omp_set_dynamic(0);
    omp_set_num_threads(10);
    #pragma omp parallel for
    for(ii = 0; ii < (int)parameter_vec(3); ii++){

    	if(ii % 50 == 0){
    		cout << "DisorderLoop at i = " << ii << endl;
    	}

    	//-----Building 1-Exciton Matrices-----
    	Matrix<double> evecs;
    	Vektor<double> evals;

    	Vektor<double> sitevec_random;
    	sitevec_random = sitevec;

    	//-----The boost Normaldistribution has the inputform (MEAN-Value, Sigma)-----
    	//-----parameter_vec(2)/2.35482 contains FWHM = 170cm^-1-----
		#pragma omp critical
        {
    		RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
        }
    	addsiteenergies(excmatrix_domains, domains, sitevec_random);
    	diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(0));

    	//-----Building 2-Exciton Matrices-----
    	Matrix<double> Exc2Matrix;
    	Matrix<double> evecs2N;
    	Vektor<double> evals2N;
    	Matrix<double> Exc2_Map;

    	Build_2Exiton_Matrix_Shifts(excmatrix_domains(0), Exc2Matrix, Exc2Shifts);
//    	Build_2Exiton_Matrix_DEBUG(excmatrix_domains(0), Exc2Matrix, Exc2Shifts);

    	diagonalizeexcmatrix_domains(evecs2N, evals2N, chlname_vec, parameter_vec, Exc2Matrix);

    	Build_2Exc_Evec_Map(excmatrix_domains(0), Exc2_Map);

//    	cout << Exc2_Map << endl;
//    	cout << "Evals2Exc: " << evals2N << "   Evecs(0,0): " << evecs2N(0,0) << endl;
//    	cout << "Evals2Exc size(): " << evals2N.size() << " evecsizes: " << evecs2N.size1() << "  " << evecs2N.size2() << endl;
//    	cout << evals2N(Exc2_Map(0,1)) << " " << evals2N(Exc2_Map(0,2)) << " " << evals2N(Exc2_Map(1,2)) << endl;
    	//double y_2N2K (int N2, int K2, Matrix<double> evec2Exc, Matrix<double> Exc2_Map, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat){
    	//
    	//    cout << "y2N2N: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "w2N0 : " << w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat) << endl;
    	//
    	//    double w2N0 = w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat);
    	//
    	//    cout << "1-Exciton Evecs: " << evecs << endl;
    	//
    	//    cout << "y2N2NMM: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "y2N2NMM: " << y_2N2KMM(0, 1, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//
    	//    cout << "tau2N: " << tau2N_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//
    	//    cout << "G(t) Prefactor: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
    	//    cout << "G(t) Prefactor: " << y_MNKL(0, 0, 0, 0, evecs, parameter_vec, dipolemat) << endl;
    	//    cout << "G(t) Prefactor: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;

    	//-----Initialize the Transition Dipole Moments-----
    	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
    	ublas::zero_vector<double> zero(3);

    	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {

    		trans_dipole_moment(M) = zero;
    	}

    	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {
    		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {

    			trans_dipole_moment(M) += dipolemat(0)(j) * evecs(M, j);
    		}
    	}

    	Matrix< Vektor<double> > trans_dipole_moment_2N( (int)evals2N.size(), (int)evals.size() );

    	for (int M = 0; M < (int) evals.size(); M++) {

    		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

    			trans_dipole_moment_2N(N2, M) = zero;
    		}
    	}

    	for (int M = 0; M < (int) evals.size(); M++) {
    		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

    			for (int k = 0; k < (int) evals.size(); k++) {
    				for (int l = 0; l < (int) evals.size(); l++) {

    					if (l > k) {

    						trans_dipole_moment_2N(N2, M) += evecs2N(N2, Exc2_Map(k, l)) * ( dipolemat(0)(k) * evecs(M,l) + dipolemat(0)(l) * evecs(M,k) );
    						//						cout << "trans_dipole_moment_2N" << trans_dipole_moment_2N(N2, M) << " evecs2N " << evecs2N(N2, Exc2_Map(k, l)) << " mu_k: " <<  dipolemat(0)(k) * evecs(M,k) << " mu_l: " << dipolemat(0)(l) * evecs(M,l) << endl;
    					}

    				}
    			}
    		}

    	}

    	//-----This is for the Paper; Rotational Strength-----
        Vektor<double> RM(evals.size());
        for(int i = 0; i < (int)RM.size(); i++){
            RM(i) = 0;
        }

        Vektor<double> CrossResult(3);
        Vektor<double> DMReals;

    	for(int M = 0; M < (int) RM.size(); M++){
    		for(int i = 0; i < (int) evals.size() - 1; i++){
    			for(int j = i + 1; j < (int) evals.size(); j++){

                    crossprod(dipolemat(0)(i), dipolemat(0)(j),CrossResult);
                    RM(M) += 1.0 * evecs(M,i) * evecs(M,j) * inner_prod( dipolemat(1)(i) - dipolemat(1)(j) , CrossResult);
                }
            }
        }

    	//-----Initial Population PM(0)-----
    	Vektor<double> P0;
    	Vektor<double> U_K;

    	//-----Simple Dimer-----
//    	push_back(P0, 1.0);
//    	push_back(P0, 0.0);
//    	cout << "evals & p0: " << endl;
//    	cout << evals << endl;
//    	cout << P0 << endl;

    	double U_w = 0.0;
    	Vektor<double> DMRealsP0;
    	fftw_complex *in3, *out3;

    	in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
    	out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

    	for (int M = 0; M < (int) evals.size(); M++) {

    		U_w = PM0_WZ(M, DMRealsP0, GReals, GImags, 8192.0 / 2048.0, ntime, CIm, evecs, evals, parameter_vec, dipolemat, in3, out3, parameter_vec(6), parameter_vec(5));
    		//push_back(P0, U_w * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) );
    		push_back(U_K, U_w);
    		push_back(P0, U_w * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) );
    	}

    	fftw_free(in3);
    	fftw_free(out3);

    	U_K = U_K / norm_2(U_K);
    	P0 = P0 / norm_2(P0);

//    	//-----DEBUG-----
//    	cout <<"DEBUG -- P0: " << P0 << " " << P0 / norm_2(P0) << endl;
////    	cout <<"DEBUG -- U_K: " << U_K << endl;
////    	cout << "mu_0: " << inner_prod(trans_dipole_moment(0), trans_dipole_moment(0)) << endl;
////    	cout << "mu_1: " << inner_prod(trans_dipole_moment(1), trans_dipole_moment(1)) << endl;
////    	exit(1);
//
//    	//	//-----DEBUG-----
//    	cout <<"DEBUG -- Evals / Evecs: " << evals << endl;
//    	//	//-----DEBUG-----
//    	cout <<"DEBUG -- Evals2N / Evecs2N: " << evals2N << endl << evecs2N << endl;
//    	//
    	//	//-----DEBUG-----

    	for(int M = 0; M < (int)trans_dipole_moment.size(); M++){
    		cout << "M : " << M << " " << trans_dipole_moment(M) << " " << norm_2(trans_dipole_moment(M)) << endl;
    	}

    	for(int M = 0; M <(int)RM.size(); M++){
    		cout << "M : " << M << " " << RM(M) << endl;
    	}

    	exit(1);
//    	cout <<"DEBUG -- Transition Dipoles: " << trans_dipole_moment << endl;
    	//
    	//	//-----DEBUG-----
    	//	cout <<"DEBUG -- Transition Dipoles 2N: " << trans_dipole_moment_2N << endl;

    	//-----Solve Rate Equation for Kinetic Matrix A-----

    	//-----Build Kinetic Matrix A-----
    	Matrix<double> A_kinetic((long)sitevec.size(), (long)sitevec.size());
    	for (long M = 0; M < (long) A_kinetic.size1(); M++) {
    		for (long N = 0; N < (long) A_kinetic.size2(); N++) {

    			A_kinetic(M, N) = 0;
    		}
    	}

		for(int M = 0; M < (int)evals.size(); M++){
			for(int N = 0; N < (int)evals.size(); N++){

				if(M != N){
					double rate = -2.0 * y_MN(M, N, evecs) * C_Re_WZ(evals(M) - evals(N), parameter_vec(0) );
					A_kinetic(M, N) = rate;
				}
			}
		}

		auto A_kinetic_fin = trans(A_kinetic);
//		cout << A_kinetic_fin << endl;

		//-----Sum Up the Diagonals-----

//		auto A_kinetic_fin = A_kinetic;
		for(long N = 0; N < (long)A_kinetic_fin.size1(); N++){
			double k_sum = 0.0;
			for(long M = 0; M < (long)A_kinetic_fin.size2(); M++){

					if(M != N){
						k_sum += A_kinetic_fin(M, N);
					}
				}
				A_kinetic_fin(N, N) = -k_sum;
		}

		A_kinetic = -A_kinetic_fin;
//		cout << "A_kinetic : " << A_kinetic << endl;

		//-----Equilibrium Populations-----
		Vektor<double> P_eq(A_kinetic.size1());

		auto P_boltzman = [&](double wM){
			return exp((-6.582*pow (10.0, -16)*1.88496*pow (10.0, 11) * wM)/(8.617*pow (10, -5)*parameter_vec(0)));
		};

		long double P_boltzman_sum = 0;

		for(long j = 0; j < (long)evals.size(); j++){
			P_boltzman_sum += P_boltzman(evals(j));
		}

		for(long j = 0; j < (long)evals.size(); j++){
			P_eq(j) = P_boltzman(evals(j)) / P_boltzman_sum;
		}

		double P_eq_sum = 0;
		for(long j = 0; j < (long)P_eq.size(); j++){
			P_eq_sum += P_eq(j);
		}

//		cout << "P_eq: " << P_eq << " P_eq_sum: " << P_eq_sum << endl;

		//-----Build Symmetic A_kinetic-----
		auto A_kinetic_sym = A_kinetic;

		for(long M = 0; M < (long)A_kinetic.size1(); M++){
			for(long N = 0; N < (long)A_kinetic.size2(); N++){

				A_kinetic_sym(M, N) = A_kinetic(M, N) * P_eq(N) / sqrt(P_eq(M) * P_eq(N));
			}
		}

//		cout << A_kinetic_sym << endl;

		Matrix<double> evecs_kinetic;
		Vektor<double> evals_kinetic;

//		for(long i = 0; i < (long) evecs_kinetic.size1(); i++){
//			cout << "{";
//			for(long j = 0; j < (long) evecs_kinetic.size2(); j++){
//				cout << evecs_kinetic(i, j) << ",";
//			}
//			cout << "}" << endl;
//		}

		diagonalize_KineticMatrix_sym(A_kinetic_sym, evecs_kinetic, evals_kinetic);
//		diagonalize_KineticMatrix(A_kinetic, evecs_kinetic, evals_kinetic);

		//-----Set Up the C_i's -> it's a Vektor of Vektors-----
		Vektor< Vektor<double> > ci_evecs;
		ci_evecs.resize((long)evecs_kinetic.size1());

		//-----Resize & Copy-----
		for(long i = 0; i < (long) ci_evecs.size(); i++){
			ci_evecs(i).resize((long)evecs_kinetic.size2());
		}

		for(long i = 0; i < (long) evecs_kinetic.size1(); i++){
			for(long j = 0; j < (long) evecs_kinetic.size2(); j++){
				ci_evecs(i)(j) = evecs_kinetic(i, j);
			}
		}
		//-----Normalize to 1 & Print-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){

			ci_evecs(i) = ci_evecs(i) / norm_2(ci_evecs(i));
//			cout << ci_evecs(i) / norm_2(ci_evecs(i)) << endl;
		}
		//-----Copy back the normalized Vektors to initial evecs_kinetic Matrix-----
		for(long i = 0; i < (long)ci_evecs.size(); i++){
			for(long j = 0; j < (long)ci_evecs.size(); j++){

				evecs_kinetic(i, j) = ci_evecs(i)(j);
			}
		}

		//-----Delay Time conversion fs -> cm-----
		double delay = parameter_vec(7) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);
//		delay = 1;
		Vektor<double> PM_t(evals_kinetic.size());

		for(int M = 0; M < (int)PM_t.size(); M++){
			PM_t(M) = 0;

			for(int i = 0; i < (int)ci_evecs.size(); i++){
				for(int N = 0; N < (int)P_eq.size(); N++){

					PM_t(M) += P0(N) / P_eq(N) * ci_evecs(i)(N) * sqrt(P_eq(N)) * ci_evecs(i)(M) * sqrt(P_eq(M)) * exp(delay * evals_kinetic(i));
				}
			}
		}

//		//-----This is for the Paper-----
//
//		ofstream output_Pt;
//		output_Pt.open("Pt_hom.out");
//
//		Vektor<double> PM_t_out(evals_kinetic.size());
////		for(int M = 0; M < (int)PM_t_out.size(); M++){
////			PM_t_out(M) = 0;
////		}
//
//		for(double delay_t = 0; delay_t < 5; delay_t = delay_t + 0.0001){
//
//			output_Pt << delay_t / ((3.0 * 1e10) * 1e-15 * (2 * 3.1415926)) << " ";
//
//			for(int M = 0; M < (int)PM_t.size(); M++){
//
//				PM_t_out(M) = 0;
//
//				for(int i = 0; i < (int)ci_evecs.size(); i++){
//					for(int N = 0; N < (int)P_eq.size(); N++){
//
//						PM_t_out(M) += P0(N) / P_eq(N) * ci_evecs(i)(N) * sqrt(P_eq(N)) * ci_evecs(i)(M) * sqrt(P_eq(M)) * exp(delay_t * evals_kinetic(i));
//					}
//				}
//
//				output_Pt << PM_t_out(M) << " ";
//			}
//			output_Pt << endl;
//		}
//
//		output_Pt.close();
//		cout << P0 << endl;
//		cout << P_eq << endl;
//
//		exit(1);

//		cout << "evals: " << evals << endl;
//		cout << P0 << endl;
//		cout << delay << endl;
//		cout << PM_t << endl;

//		//-----Check that W_KM works-----
//		Vektor<double> PM_t_WK(evals_kinetic.size());
//
//		for(int M = 0; M < (int)PM_t.size(); M++){
//			PM_t_WK(M) = 0;
//			for(int K = 0; K < (int) PM_t.size(); K++){
//				PM_t_WK(M) += W_KM(K, M, delay, ci_evecs, evals_kinetic, P_eq, P0);
//			}
//		}
//
//		cout << PM_t_WK << endl;

//		cout << "Normierte Kinetic Matrix: " << evecs_kinetic << endl;
//		cout << "evals_kineitc: " << evals_kinetic << endl;
//
//		//-----Get Initial conditions d_i from P(0)-----
//		//-----Solve the system of linear equations-----
//		Matrix<double> c_matrix;
//		c_matrix = trans(evecs_kinetic);
//
//		Vektor<double> d;
//		d = solve_lineq(c_matrix, P0);
//
//		cout << "c_matrix: " << c_matrix << endl;
//		cout << "Initidal Conditions d_i: " << d << endl;
//		cout << "starting from Vektor P(0) = " << P0 << endl;
//		cout << "evals_kinetic: " << evals_kinetic << endl;
//		cout << "evecs_kinetic: " << evecs_kinetic << endl;
		//-----Output P(t)-----
//		double sum = 0;
//
//		ofstream Pt;
//		Pt.open("PT_Vergleich_Problem.out");
//
//		double t_max = 150;
//		double dt = 0.005;
//		int num_steps = (int)(t_max / dt);
//		Vektor< Vektor<double> > population(num_steps);
//
//		for(int k = 0; k < (int)population.size(); k++){
//			population(k).resize(evals.size());
//		}
//
//		for(long k = 0; k < num_steps; k++){
//			Pt << k * dt << " ";
//			for(int i = 0; i < (int)d.size(); i++){
//
//				population(k) += d(i) * ci_evecs(i) * exp(dt * k * evals_kinetic(i));
//			}
//			Pt << population(k) << endl;
//		}
//
//		Pt.close();
//
//		//-----Delay Time conversion fs -> cm-----
//		double delay = parameter_vec(7) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);
//		delay = 0;
//		Vektor<double> PM_t(evals_kinetic.size());
//
//		for(int i = 0; i < (int)d.size(); i++){
//			PM_t += d(i) * ci_evecs(i) * exp(delay * evals_kinetic(i));
//		}
//
//		cout << delay << endl;
//		cout << PM_t << endl;
//
//		//-----Equilibrium Populations-----
//
//		double t_inf = 50000 * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);
//		Vektor<double> PM_equ(evals_kinetic.size());
//		Vektor<double> PM_equ_2(evals_kinetic.size());
//
//		for(int i = 0; i < (int)d.size(); i++){
//			PM_equ += d(i) * ci_evecs(i) * exp(t_inf * evals_kinetic(i));
//		}
//
//		for(int i = 0; i < (int)d.size(); i++){
//			PM_equ_2 += d(i) * ci_evecs(i) * exp(t_inf * 2 * evals_kinetic(i));
//		}
//
//		//-----Check if equilibrium is reached-----
//		double error = 0;
//		for(int M = 0; M < (int)PM_equ.size(); M++){
//			error += abs(PM_equ(M) - PM_equ_2(M));
//		}
//
//		if(error > 1e-3){
//			cout << "DEBUG: Equlibration-Error: " << error << endl;
//			cout << "Check Equilibration-Time!! Exiting for now!" << endl;
//			exit(1);
//		}
//
//		cout << "DEBUG: Equilibrium Population: " << PM_equ << endl;
//
//		cout << "PM(delay) = " << PM_t << endl;
//
//		Vektor<double> PM_delay_WKM(evals_kinetic.size());
//
//		for(int M = 0; M < (int)evals_kinetic.size(); M++){
//			PM_delay_WKM(M) = 0;
//			for(int j = 0; j < (int)ci_evecs.size(); j++){
//				for(int N = 0; N < (int)evals_kinetic.size(); N++){
//	//				PM_delay_WKM(M) += inner_prod(trans_dipole_moment(K), trans_dipole_moment(K)) * W_KM(K, M, delay, ci_evecs, evals_kinetic, PM_equ, U_vec);
//
//					PM_delay_WKM(M) +=	U_vec(N) / PM_equ(N) * inner_prod(trans_dipole_moment(M), trans_dipole_moment(M)) * ci_evecs(j)(N) * ci_evecs(j)(M) * exp(delay * evals_kinetic(j));
//				}
//			}
//		}
//
//		double PM_delay_symm = 0;
//
//		for(int j = 0; j < (int)ci_evecs.size(); j++){
//
//				PM_delay_symm += 1.0 / PM_equ(0) * ci_evecs(j)(0) * ci_evecs(j)(0) * exp(evals_kinetic(j));
//				cout << PM_delay_symm << " " << 1.0 / PM_equ(0) << " " << ci_evecs(j)(0) * ci_evecs(j)(0) << " " << exp(evals_kinetic(j)) << endl;
//		}
//
//		cout << "PM_WKM(delay) = " << PM_delay_symm << endl;
//		cout << evals_kinetic << endl;

//		//-----Time-Resolved Part-----
//		//-----Groundstate Bleaching-----
//    	Vektor<double> GB;
//    	fftw_complex *inGB, *outGB;
////    	inGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
////    	outGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//
//    	Lin_GSBleach(evecs, evals, dipolemat, parameter_vec, GB, CIm, GReals, GImags, inGB, outGB, P0, trans_dipole_moment);
//
////    	//-----GB Output Stuff-----
////    	FILE * ODdisorder;
////    	string ODFileName = "GB_Test.out";
////    	ODdisorder = fopen(ODFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(ODdisorder, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , GB(w));
////    	}
////    	fclose(ODdisorder);
//
//    	//-----Stimulated Emission-----
//    	Vektor<double> SE;
//    	fftw_complex *inSE, *outSE;
////    	inSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
////    	outSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//
//    	Lin_StimEmission(evecs, evals, dipolemat, parameter_vec, SE, CIm, GReals, GImags, inSE, outSE, P0, trans_dipole_moment, PM_t);
////    	Lin_StimEmission_WKM(evecs, evals, dipolemat, parameter_vec, SE, CIm, GReals, GImags, inSE, outSE, trans_dipole_moment, P0, P_eq, evals_kinetic, ci_evecs, delay);
//
////    	//-----SE Output Stuff-----
////    	FILE * SEDis;
////    	string SEFileName = "SE_Test.out";
////    	SEDis = fopen(SEFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(SEDis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , SE(w));
////    	}
//
//    	//-----Excited State Absorption-----
//    	Vektor<double> ESA;
//    	fftw_complex *inESA, *outESA;
////    	inESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
////    	outESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//
//    	Lin_ESA(evecs, evals, evecs2N, evals2N, Exc2_Map, dipolemat, parameter_vec, ESA, CIm, GReals, GImags, inESA, outESA, P0, trans_dipole_moment, trans_dipole_moment_2N, PM_t);
//
////    	//-----ESA Output Stuff-----
////    	FILE * ESADis;
////    	string ESAFileName = "ESA_Test.out";
////    	ESADis = fopen(ESAFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(ESADis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , ESA(w));
////    	}
//
//    	//-----Linear Part-----
    	//------Linear OD & CD Spectra-----
    	Vektor<double> CD;
    	Vektor<double> alpha;

    	CD_OD_Markov(evecs, evals, dipolemat, parameter_vec, CD, alpha, CIm, trans_dipole_moment);

//    	//-----CD & OD Output Stuff-----
//    	FILE * CDODDis;
//    	string CDODFileName = "CD_OD_Test.out";
//    	CDODDis = fopen(CDODFileName.c_str(), "wb");
//
//    	for (int w = 0; w < 15000; w++) {
//    		fprintf(CDODDis, " %f  %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , CD(w), alpha(w));
//    	}
//
//    	fclose(CDODDis);

//    	//-----TRCD Spectra-----
//    	//-----CD Groundstate Bleaching-----
//    	Vektor<double> GBCD;
//    	CD_GSBleach(evecs, evals, dipolemat, parameter_vec, GBCD, CIm, trans_dipole_moment, P0, P_eq, evals_kinetic, ci_evecs, delay);
//
////    	//-----GB Output Stuff-----
////    	FILE * CDGSBdis;
////    	string CDGSBFileName = "CDGB_Test.out";
////    	CDGSBdis = fopen(CDGSBFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(CDGSBdis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , GBCD(w));
////    	}
////    	fclose(CDGSBdis);
////
//    	//-----CD Stimulated Emission-----
//    	Vektor<double> SECD;
//    	CD_SE(evecs, evals, dipolemat, parameter_vec, SECD, CIm, trans_dipole_moment, P0, P_eq, evals_kinetic, ci_evecs, delay);
//
////    	//-----SE Output Stuff-----
////    	FILE * CDSEdis;
////    	string CDGSEFileName = "CDSE_Test.out";
////    	CDSEdis = fopen(CDGSEFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(CDSEdis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , SECD(w));
////    	}
////    	fclose(CDSEdis);
////
//    	//-----CD Excited State Absorption-----
//    	Vektor<double> ESACD;
//    	CD_ESA(evecs, evals, evecs2N, evals2N, Exc2_Map, dipolemat, parameter_vec, ESACD, CIm, trans_dipole_moment, P0, P_eq, evals_kinetic, ci_evecs, delay);
//
////    	//-----ESA Output Stuff-----
////    	FILE * CDESAdis;
////    	string CDESAFileName = "CDESA_Test.out";
////    	CDESAdis = fopen(CDESAFileName.c_str(), "wb");
////
////    	for (int w = 0; w < 15000; w++) {
////    		fprintf(CDESAdis, " %f  %f\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , ESACD(w));
////    	}
////    	fclose(CDESAdis);
//
//    	//----part2----

    	#pragma omp critical
    	{
//    		GB_sum(ii) = GB;
//    		SE_sum(ii) = SE;
//    		ESA_sum(ii) = ESA;
    		CD_sum(ii) = CD;
    		OD_sum(ii) = alpha;
//    		CDGB_sum(ii) = GBCD;
//    		CDSE_sum(ii) = SECD;
//    		CDESA_sum(ii) = ESACD;
    	}

    	sitevec_random = sitevec;

//    	fftw_free(inGB);
//    	fftw_free(outGB);
//
//    	fftw_free(inSE);
//    	fftw_free(outSE);
//
//    	fftw_free(inESA);
//    	fftw_free(outESA);
    }

    Vektor<double> GB_sum_total;
	Vektor<double> SE_sum_total;
	Vektor<double> ESA_sum_total;
	Vektor<double> Pump_Test_total;

	Vektor<double> CD_sum_total;
	Vektor<double> OD_sum_total;

    Vektor<double> CDGB_sum_total;
	Vektor<double> CDSE_sum_total;
	Vektor<double> CDESA_sum_total;
	Vektor<double> TRCD_total;

	GB_sum_total.resize(GB_sum(0).size());
	SE_sum_total.resize(SE_sum(0).size());
	ESA_sum_total.resize(ESA_sum(0).size());
	Pump_Test_total.resize(ESA_sum(0).size());

	CD_sum_total.resize(CD_sum(0).size());
	OD_sum_total.resize(OD_sum(0).size());

	CDGB_sum_total.resize(CDGB_sum(0).size());
	CDSE_sum_total.resize(CDSE_sum(0).size());
	CDESA_sum_total.resize(CDESA_sum(0).size());
	TRCD_total.resize(CDESA_sum(0).size());

	for (int j = 0; j < (int) GB_sum.size(); j++) {
		GB_sum_total = GB_sum_total + GB_sum(j);
		SE_sum_total = SE_sum_total + SE_sum(j);
		ESA_sum_total = ESA_sum_total + ESA_sum(j);
		CD_sum_total = CD_sum_total + CD_sum(j);
		OD_sum_total = OD_sum_total + OD_sum(j);
		CDGB_sum_total = CDGB_sum_total + CDGB_sum(j);
		CDSE_sum_total = CDSE_sum_total + CDSE_sum(j);
		CDESA_sum_total = CDESA_sum_total + CDESA_sum(j);
	}

	GB_sum_total = GB_sum_total / parameter_vec(3);
	SE_sum_total = SE_sum_total / parameter_vec(3);
	ESA_sum_total = ESA_sum_total / parameter_vec(3);
	CD_sum_total = CD_sum_total / parameter_vec(3);
	OD_sum_total = OD_sum_total / parameter_vec(3);
	CDGB_sum_total = CDGB_sum_total / parameter_vec(3);
	CDSE_sum_total = CDSE_sum_total / parameter_vec(3);
	CDESA_sum_total = CDESA_sum_total / parameter_vec(3);

	Pump_Test_total = GB_sum_total + SE_sum_total + ESA_sum_total;
	TRCD_total = CDGB_sum_total + CDSE_sum_total + CDESA_sum_total;

	cout << "Calculations done, writing out files!" << endl;

//	FILE * Pump_Test_Disorder;
////	string Pump_Test_FileName = "Pump_Test_Dis.out";
//	string Pump_Test_FileName = "Pump_Test_Dis_" + to_string((int)parameter_vec(7)) + ".out";
//
//	Pump_Test_Disorder = fopen(Pump_Test_FileName.c_str(), "wb");
//
//	for (int w = 0; w < 15000; w++) {
//		fprintf(Pump_Test_Disorder, " %lf  %lf  %lf  %lf  %lf\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , Pump_Test_total(w), GB_sum_total(w), SE_sum_total(w), ESA_sum_total(w));
//	}
//	fclose(Pump_Test_Disorder);
//
//	FILE * TRCD_Disorder;
//	string TRCD_FileName = "TRCD_Dis_" + to_string((int)parameter_vec(7)) + ".out";
//
//	TRCD_Disorder = fopen(TRCD_FileName.c_str(), "wb");
//
//	for (int w = 0; w < 15000; w++) {
//		fprintf(TRCD_Disorder, " %lf  %lf  %lf  %lf  %lf\n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , TRCD_total(w), CDGB_sum_total(w), CDSE_sum_total(w), CDESA_sum_total(w));
//	}
//	fclose(TRCD_Disorder);

	FILE * CD_OD_Disorder;
	string CD_OD_FileName = "CD_OD_Dis.out";

	CD_OD_Disorder = fopen(CD_OD_FileName.c_str(), "wb");
	for (int w = 0; w < 15000; w++) {
		fprintf(CD_OD_Disorder, " %lf  %lf  %lf \n", 1.0 / (1.0/(1000 * 1e-7) + w) * 1e7 , CD_sum_total(w), OD_sum_total(w));
	}
	fclose(CD_OD_Disorder);

	cout << "Writing files done!" << endl;
}

//void Lin_Spectra_Markov(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, double t_ang){
//
//	//-----Linear OD and CD with Static Disorder-----
//	//-----Initializations-----
//
//    Vektor< Vektor<double> > alpha_sum ((int)parameter_vec(3));
//    Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
//    Vektor< Vektor<double> > RotStrength_sum ((int)parameter_vec(3));
//
//    Vektor< Matrix<double> > excmatrix_domains;
//    Vektor< Vektor<double> > domains;
//
//    //Ausgliederung von CIm & G(t)
//    Vektor<double> CIm;
//    correlation(parameter_vec, CIm);
//
//    Vektor<double> GReals;
//    Vektor<double> GImags;
//    FT_Gt_WZ(GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,parameter_vec);
//
//    Vektor< Vektor<double> > coupling_vec;
//    buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);
//
//    double ntime = 8192.0 * 4.0;
//
//    printf("Initialization Done! Doing Static-Disorder Loop now.\n");
//
//    int i = 0;
////	#pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
//	#pragma omp parallel for
//    for(i = 0; i < (int)parameter_vec(3); i++){
//
//        Vektor<double> CD;
//        Vektor<double> alpha;
//        Matrix<double> evecs;
//        Vektor<double> evals;
//
//        Vektor<double> RotStrength_vec;
//
//        Vektor<double> sitevec_random;
//        sitevec_random = sitevec;
//
//        fftw_complex *in, *out;
//        fftw_plan p;
//
//        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
//        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
//
//        RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
////        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
//        addsiteenergies(excmatrix_domains,domains,sitevec_random);
//
//        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;
//
//        for(int j = 0; j < (int)excmatrix_domains.size(); j++){  //(int)excmatrix_domains.size()
//
//            Vektor< Vektor< Vektor<double> > > dipolemat_domains;
//            dipolemat_domains = dipolemat;
//
//            dipolemat_domains(0).resize(domains(j).size());
//            dipolemat_domains(1).resize(domains(j).size());
//
//            for(int u = 0; u < (int)dipolemat_domains(0).size(); u++){
//
//                dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u)-1);
//                dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u)-1);
//            }
//
//            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));
////            get_RotStrength(evecs, evals, dipolemat_domains, parameter_vec, RotStrength_vec, 14814.8);
//
//            CDabsorption_domains(evecs, evals, dipolemat_domains, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);
//
//
//			#pragma omp critical
//            {
//                CD_sum(i) = CD;
//                alpha_sum(i) = alpha;
//                RotStrength_sum(i) = RotStrength_vec;
//            }
//            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
//        }
//
//        sitevec_random = sitevec;
//        fftw_destroy_plan(p);
//        fftw_free(in);
//        fftw_free(out);
//    }
//
////    FILE * RotStrength_dis_out;
////    RotStrength_dis_out = fopen ("RotStrength.out", "a+");
////
////    fprintf(RotStrength_dis_out,"<R->dis = %lf at Angle %lf \n", dis_RotStrength(RotStrength_sum, -1), t_ang);
////    cout << "dis_RotStrength: " << dis_RotStrength(RotStrength_sum, -1) << endl;
////
////    fclose(RotStrength_dis_out);
//
//    Vektor<double> CD_sum_total;
//    Vektor<double> alpha_sum_total;
//
//    CD_sum_total.resize(CD_sum(0).size());
//    alpha_sum_total.resize(alpha_sum(0).size());
//
//    for(int i = 0; i < (int)CD_sum.size(); i++){
//        CD_sum_total = CD_sum_total + CD_sum(i);
//        alpha_sum_total = alpha_sum_total + alpha_sum(i);
//    }
//
//    alpha_sum_total = alpha_sum_total / parameter_vec(3);
//    CD_sum_total = CD_sum_total / parameter_vec(3);
//
//    cout<<"Calculations done, writing out files!"<<endl;
//
//    FILE * ODdisorder;
//    FILE * CDdisorder;
//
//    string ODFileName = "OD_Dis.out";
//    string CDFileName = "CD_Dis.out";
//
//    ODdisorder = fopen ( ODFileName.c_str() , "wb");
//    CDdisorder = fopen ( CDFileName.c_str() , "wb");
//
//    for(int u = 0; u < (int)CD_sum_total.size(); u++){
//
//        fprintf(ODdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),alpha_sum_total(u));  //16
//        fprintf(CDdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),CD_sum_total(u));     //16
//    }
//
//    fclose (ODdisorder);
//    fclose (CDdisorder);
//
//    cout<<"Writing files done!"<<endl;
//
//    return 0;
//}

void CD_Pump_Test_Dimerold(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat){
	//-----CD Pump-Test with Static Disorder-----

	//-----Initializations-----

	Vektor<Vektor<double> > GB_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > SE_sum((int) parameter_vec(3));
	Vektor<Vektor<double> > ESA_sum((int) parameter_vec(3));

	double ntime = 8192.0 * 4.0;

	//Ausgliederung von CIm & G(t)

	Vektor<double> GReals;
	Vektor<double> GImags;
	FT_Gt_WZ(GReals, GImags, 8192.0 / 2048.0, ntime, parameter_vec); //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0

	//cout<<"Domainentry :  "<<domains(0)(0)<<endl;
	//cout<<endl<<endl<<"Dipolemat :  "<<dipolemat<<endl;
	//cout<<"Excmatrix_domains size :  "<<excmatrix_domains.size()<<endl;

	Vektor<Matrix<double> > excmatrix_domains;
	Vektor<Vektor<double> > domains;

	Vektor<double> CIm;						//Ausgliederung von CIm & G(t)
	correlation(parameter_vec, CIm);

	Vektor<Vektor<double> > coupling_vec;
	buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    printf("Initialization for Pump-Test Done! Doing Static-Disorder Loop now.\n");

    //-----Insert Static-Disorder-----

    int i = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
    for(i = 0; i < (int)parameter_vec(3); i++){ //parameter_vec(3) holds #nrand


    if(i % 40 == 0){
    	cout << "DisorderLoop at i = " << i << endl;
    }

	//-----Building 1-Exciton Matrices-----
    Matrix<double> evecs;
    Vektor<double> evals;

    Vektor<double> sitevec_random;
    sitevec_random = sitevec;

    //-----The boost Normaldistribution has the inputform (MEAN-Value, Sigma)-----
    //-----parameter_vec(2)/2.35482 contains FWHM = 170cm^-1-----

    RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
    addsiteenergies(excmatrix_domains, domains, sitevec_random);

    diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(0));

	//-----Building 2-Exciton Matrices-----
    Matrix<double> Exc2Matrix;
    Matrix<double> evecs2N;
    Vektor<double> evals2N;
    Matrix<double> Exc2_Map;

    Build_2Exiton_Matrix(excmatrix_domains(0), Exc2Matrix);

    diagonalizeexcmatrix_domains(evecs2N, evals2N, chlname_vec, parameter_vec, Exc2Matrix);

    Build_2Exc_Evec_Map(excmatrix_domains(0), Exc2_Map);

//    cout <<"Exc2 Matrix : " << Exc2Matrix << endl;
//    cout << "Evals2Exc: " << evals2N << "   Evecs(0,0): " << evecs2N(0,0) << endl;
//    cout << "Evals2Exc size(): " << evals2N.size() << " evecsizes: " << evecs2N.size1() << "  " << evecs2N.size2() << endl;
//    //double y_2N2K (int N2, int K2, Matrix<double> evec2Exc, Matrix<double> Exc2_Map, Vektor<double> parameter_vec, Vektor< Vektor< Vektor<double> > > dipolemat){
//
//    cout << "y2N2N: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
//    cout << "w2N0 : " << w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat) << endl;
//
//    double w2N0 = w2N0_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, CIm, dipolemat);
//
//    cout << "1-Exciton Evecs: " << evecs << endl;
//
//    cout << "y2N2NMM: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
//    cout << "y2N2NMM: " << y_2N2KMM(0, 1, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;
//
//    cout << "tau2N: " << tau2N_WZ(0, evals2N, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
//
//    cout << "G(t) Prefactor: " << y_2N2K(0, 0, evecs2N, Exc2_Map, parameter_vec, dipolemat) << endl;
//    cout << "G(t) Prefactor: " << y_MNKL(0, 0, 0, 0, evecs, parameter_vec, dipolemat) << endl;
//    cout << "G(t) Prefactor: " << y_2N2KMM(0, 0, evecs2N, evecs, Exc2_Map, parameter_vec, dipolemat) << endl;


	//-----Initialize the Transition Dipole Moments-----
	Vektor<Vektor<double> > trans_dipole_moment(evals.size());
	ublas::zero_vector<double> zero(3);

	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {

		trans_dipole_moment(M) = zero;
	}

	for (int M = 0; M < (int) trans_dipole_moment.size(); M++) {
		for (int j = 0; j < (int) trans_dipole_moment.size(); j++) {

			trans_dipole_moment(M) += dipolemat(0)(j) * evecs(M, j);      //dipolemat(0) because it contains the ei_dipole
		}
	}

	Matrix< Vektor<double> > trans_dipole_moment_2N( (int)evals2N.size(), (int)evals.size() );

	for (int M = 0; M < (int) evals.size(); M++) {

		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

			trans_dipole_moment_2N(N2, M) = zero;
		}
	}


	for (int M = 0; M < (int) evals.size(); M++) {
		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

			for (int k = 0; k < (int) evals.size(); k++) {
				for (int l = 0; l < (int) evals.size(); l++) {

					if (l > k) {

						trans_dipole_moment_2N(N2, M) += evecs2N(N2, Exc2_Map(k, l)) * ( dipolemat(0)(k) * evecs(M,l) + dipolemat(0)(l) * evecs(M,k) );
//						cout << "trans_dipole_moment_2N" << trans_dipole_moment_2N(N2, M) << " evecs2N " << evecs2N(N2, Exc2_Map(k, l)) << " mu_k: " <<  dipolemat(0)(k) * evecs(M,k) << " mu_l: " << dipolemat(0)(l) * evecs(M,l) << endl;
					}

				}
			}
		}

	}

	//-----Magnetic Transition Dipole Moments-----
	Vektor<Vektor<double> > magn_dipole_moment(evals.size());
	ublas::zero_vector<double> zero_magn(3);
    Vektor<double> CrossResult(3);

	for (int M = 0; M < (int) magn_dipole_moment.size(); M++) {

		magn_dipole_moment(M) = zero_magn;
	}

	for (int M = 0; M < (int) magn_dipole_moment.size(); M++) {
		for (int j = 0; j < (int) magn_dipole_moment.size(); j++) {


			crossprod(dipolemat(1)(j), dipolemat(0)(j), CrossResult);
			magn_dipole_moment(M) += CrossResult * evecs(M, j);

			//	//-----DEBUG-----
//			cout << "M : " << M << "  j: " << j << "  CR: " << CrossResult << "  evecs: " << evecs(M,j) << "   m: " << magn_dipole_moment(M) << endl;
		}
	}

	//	//-----DEBUG-----
//	cout <<"DEBUG: Magnetic Moments mM " << magn_dipole_moment << endl;

	Matrix< Vektor<double> > magn_dipole_moment_2N( (int)evals2N.size(), (int)evals.size() );

	for (int M = 0; M < (int) evals.size(); M++) {

		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

			magn_dipole_moment_2N(N2, M) = zero;
		}
	}

    Vektor<double> CrossResult_k(3);
    Vektor<double> CrossResult_l(3);

	for (int M = 0; M < (int) evals.size(); M++) {
		for (int N2 = 0; N2 < (int) evals2N.size(); N2++) {

			for (int k = 0; k < (int) evals.size(); k++) {
				for (int l = 0; l < (int) evals.size(); l++) {

					if (l > k) {

						crossprod(dipolemat(1)(k) , dipolemat(0)(k), CrossResult_k);
						crossprod(dipolemat(1)(l) , dipolemat(0)(l), CrossResult_l);

						magn_dipole_moment_2N(N2, M) += evecs2N(N2, Exc2_Map(k, l)) * ( CrossResult_k * evecs(M,l) + CrossResult_l * evecs(M,k) );

//						cout << "trans_dipole_moment_2N" << trans_dipole_moment_2N(N2, M) << " evecs2N " << evecs2N(N2, Exc2_Map(k, l)) << " mu_k: " <<  dipolemat(0)(k) * evecs(M,k) << " mu_l: " << dipolemat(0)(l) * evecs(M,l) << endl;
					}

				}
			}
		}

	}

	//	//-----DEBUG-----
//	cout <<"DEBUG: 2-Exc Magnetic Moments m2N,M " << magn_dipole_moment_2N << endl;

	//-----Initial Population PM(0)-----
	Vektor<double> P0;
	double U_w = 0.0;
	Vektor<double> DMRealsP0;
	fftw_complex *in3, *out3;
	fftw_plan p3;
	in3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
	out3 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

	for (int M = 0; M < (int) evals.size(); M++) {

		U_w = PM0_WZ(M, DMRealsP0, GReals, GImags, 8192.0 / 2048.0, ntime, CIm, evecs, evals, parameter_vec, dipolemat, in3, out3, parameter_vec(6), parameter_vec(5));
		push_back(P0, U_w );
		//-----ATTENTION!: THERE THIS PART DIFFERS FROM THE LINEAR PUMP-TEST CASE AS THE |mu1 . mu1| are written explicitly into the functions below-----
	}

    fftw_destroy_plan(p3);
    fftw_free(in3);
    fftw_free(out3);


//	//-----DEBUG-----
//	cout <<"DEBUG -- P0: " << P0 << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Evals / Evecs: " << evals << endl << evecs << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Evals2N / Evecs2N: " << evals2N << endl << evecs2N << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Transition Dipoles: " << trans_dipole_moment << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Transition Dipoles 2N: " << trans_dipole_moment_2N << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Magnetic Dipoles: " << magn_dipole_moment << endl;
//
//	//-----DEBUG-----
//	cout <<"DEBUG -- Magnetic Dipoles 2N: " << magn_dipole_moment_2N << endl;


    //-----Groundstate Bleaching-----
	Vektor<double> GB;
	fftw_complex *inGB, *outGB;
	fftw_plan pGB;

	inGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
	outGB = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

	CD_GSBleach_Dimerold(evecs, evals, dipolemat, parameter_vec, GB, CIm, GReals, GImags, inGB, outGB, pGB, P0, trans_dipole_moment, magn_dipole_moment);



//	Lin_GSBleach(evecs, evals, dipolemat, parameter_vec, GB, CIm, GReals, GImags, inGB, outGB, pGB, P0, trans_dipole_moment);

//	//-----GB Output Stuff-----
//	FILE * ODdisorder;
//	string ODFileName = "GB_Test.out";
//	ODdisorder = fopen(ODFileName.c_str(), "wb");
//
//	for (int u = 0; u < (int) GB.size(); u++) {
//
//		fprintf(ODdisorder, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))),GB(u));
//	}
//	fclose(ODdisorder);


    //-----Stimulated Emission-----
	//-----Delay Time conversion fs -> cm-----
	double delay = parameter_vec(7) * (3.0 * 1e10) * 1e-15 * (2 * 3.1415926);

	Vektor<double> SE;
	fftw_complex *inSE, *outSE;
	fftw_plan pSE;

	inSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
	outSE = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

	CD_StimEmission_Dimerold(evecs, evals, dipolemat, parameter_vec, SE, CIm, GReals, GImags, inSE, outSE, pSE, P0, trans_dipole_moment, magn_dipole_moment, delay);



//	//-----SE Output Stuff-----
//	FILE * SEDis;
//	string SEFileName = "SE_Test.out";
//	SEDis = fopen(SEFileName.c_str(), "wb");
//
//	for (int u = 0; u < (int) SE.size(); u++) {
//
//		fprintf(SEDis, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), SE(u));
//	}
//	fclose(SEDis);


    //-----Excited State Absorption-----
	Vektor<double> ESA;
	fftw_complex *inESA, *outESA;
	fftw_plan pESA;

	inESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);
	outESA = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

	CD_ESA_Dimerold(evecs, evals, evecs2N, evals2N, Exc2_Map, dipolemat, parameter_vec, ESA, CIm, GReals, GImags, inESA, outESA, pESA, P0, trans_dipole_moment, trans_dipole_moment_2N, magn_dipole_moment, magn_dipole_moment_2N, delay);

//	//-----ESA Output Stuff-----
//	FILE * ESADis;
//	string ESAFileName = "ESA_Test.out";
//	ESADis = fopen(ESAFileName.c_str(), "wb");
//
//	for (int u = 0; u < (int) ESA.size(); u++) {
//
//		fprintf(ESADis, " %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ESA(u));
//	}
//	fclose(ESADis);


//	#pragma omp critical
	{
		GB_sum(i) = GB;
		SE_sum(i) = SE;
		ESA_sum(i) = ESA;
	}

    sitevec_random = sitevec;

    fftw_destroy_plan(pGB);
    fftw_free(inGB);
    fftw_free(outGB);

    fftw_destroy_plan(pSE);
    fftw_free(inSE);
    fftw_free(outSE);

    fftw_destroy_plan(pESA);
    fftw_free(inESA);
    fftw_free(outESA);

    }


    Vektor<double> GB_sum_total;
	Vektor<double> SE_sum_total;
	Vektor<double> ESA_sum_total;
	Vektor<double> Pump_Test_total;

	GB_sum_total.resize((int)GB_sum(0).size());
	SE_sum_total.resize((int)SE_sum(0).size());
	ESA_sum_total.resize((int)ESA_sum(0).size());
	Pump_Test_total.resize((int)ESA_sum(0).size());


	cout << "Sizes: " << GB_sum_total.size() << "  " << SE_sum_total.size() << "  " << ESA_sum_total.size() << "  " << Pump_Test_total.size() << endl;
	cout << "Sizes: " << GB_sum.size() << "  " << SE_sum.size() << "  " << ESA_sum.size() << "  " << endl;


	for (int j = 0; j < (int)GB_sum.size(); j++) {	//(int)GB_sum.size()
		GB_sum_total = GB_sum_total + GB_sum(j);
		SE_sum_total = SE_sum_total + SE_sum(j);
		ESA_sum_total = ESA_sum_total + ESA_sum(j);
	}

	GB_sum_total = GB_sum_total / parameter_vec(3);
	SE_sum_total = SE_sum_total / parameter_vec(3);
	ESA_sum_total = ESA_sum_total / parameter_vec(3);

	Pump_Test_total = GB_sum_total + SE_sum_total + ESA_sum_total;

	cout << "Calculations done, writing out files!" << endl;

	FILE * Pump_Test_Disorder;
	string Pump_Test_FileName = "CD_Pump_Test_Dis_" + to_string((int)parameter_vec(7)) + ".out";

	Pump_Test_Disorder = fopen(Pump_Test_FileName.c_str(), "wb");

	for (int u = 0; u < (int) Pump_Test_total.size(); u++) {

		fprintf(Pump_Test_Disorder, " %lf  %lf  %lf  %lf  %lf\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), Pump_Test_total(u), GB_sum_total(u), SE_sum_total(u), ESA_sum_total(u));  //16
		//-----w, Pump-Test, GB, SE, ESA-----
	}



	fclose(Pump_Test_Disorder);

	cout << "Writing files done!" << endl;
}

void ACD_domains_Get_Contributions(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, double t_ang, Vektor< ChlAtoms > Chl_Vec){
    //-----This Function has been used to Calculate the Baseplate ACD!-----
	//-----Static Disorder-----

    Vektor< Vektor<double> > alpha_sum_all ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_all ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Edge_sum_all ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Face_sum_all ((int)parameter_vec(3));
    Vektor< Vektor<double> > LD_sum_all ((int)parameter_vec(3));

    //-----Low Energy (LE) Contributions-----

    Vektor< Vektor<double> > alpha_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Edge_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Face_sum_LE ((int)parameter_vec(3));
    Vektor< Vektor<double> > LD_sum_LE ((int)parameter_vec(3));

    //-----High Energy (HE) Contributions-----

    Vektor< Vektor<double> > alpha_sum_HE ((int)parameter_vec(3));
    Vektor< Vektor<double> > CD_sum_HE ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Edge_sum_HE ((int)parameter_vec(3));
    Vektor< Vektor<double> > ACD_Face_sum_HE ((int)parameter_vec(3));
    Vektor< Vektor<double> > LD_sum_HE ((int)parameter_vec(3));

    Vektor< Vektor< Vektor<double> > > RotStrength_sum ((int)parameter_vec(3));
    Vektor< Vektor< Vektor<double> > > LDStrength_sum ((int)parameter_vec(3));

    Vektor< Matrix<double> > excmatrix_domains;
    Vektor< Vektor<double> > domains;

    //Ausgliederung von CIm & G(t)
    Vektor<double> CIm;
    correlation(parameter_vec, CIm);

    Vektor<double> GReals;
    Vektor<double> GImags;
    FT_Gt_WZ(GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,parameter_vec);    //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0

    Vektor< Vektor<double> > coupling_vec;
//    build_predefined_excitondomains(chlname_vec, excmatrix_domains, domains, coupling_vec);
    buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    //cout<<"Domainentry :  "<<domains(0)(0)<<endl;
    //cout<<endl<<endl<<"Dipolemat :  "<<dipolemat<<endl;
    //cout<<"Excmatrix_domains size :  "<<excmatrix_domains.size()<<endl;

    double ntime = 8192.0 * 4.0;

    printf("Initialization Done! Doing Static-Disorder Loop now.\n");

    int i = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
//    #pragma omp parallel for default(none) private(i) shared(Chl_Vec, dipolemat, ntime,excmatrix_domains, alpha_sum_all, CD_sum_all, ACD_Edge_sum_all, ACD_Face_sum_all,alpha_sum_LE, CD_sum_LE, ACD_Edge_sum_LE, ACD_Face_sum_LE, alpha_sum_HE, CD_sum_HE, ACD_Edge_sum_HE, ACD_Face_sum_HE, RotStrength_sum, parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
    for(i = 0; i < (int)parameter_vec(3); i++){

    	if(i % 200 == 0)
			cout << "Step #: " << i << endl;

        Vektor< Vektor<double> > CD(3);
        Vektor< Vektor<double> > ACD_Edge(3);
        Vektor< Vektor<double> > ACD_Face(3);
        Vektor< Vektor<double> > alpha(3);
        Vektor< Vektor<double> > LD(3);

        Matrix<double> evecs;            //matrix&vector for Eigenvectors & Eigenvalues
        Vektor<double> evals;

        Vektor< Vektor<double> > RotStrength_vec;
        Vektor< Vektor<double> > LDStrength_vec;

        Vektor<double> sitevec_random;
        sitevec_random = sitevec;

        fftw_complex *in, *out;
        fftw_plan p;

        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

        RandomizeSiteEnergies(sitevec_random, parameter_vec(2)/2.35482);
//        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
        addsiteenergies(excmatrix_domains,domains,sitevec_random);

        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;

        for(int j = 0; j < (int)excmatrix_domains.size(); j++){  //(int)excmatrix_domains.size()

            Vektor< Vektor< Vektor<double> > > dipolemat_domains;
            dipolemat_domains = dipolemat;

            //cout<<"Building dipolemat_domains :  "<<endl;

            dipolemat_domains(0).resize(domains(j).size());
            dipolemat_domains(1).resize(domains(j).size());

            for(int u = 0; u < (int)dipolemat_domains(0).size(); u++){

                dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u)-1);
                dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u)-1);
            }


            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));

//            get_RotStrength(evecs, evals, dipolemat_domains, parameter_vec, RotStrength_vec, 14814.8);

//            CDabsorption_domains(evecs, evals, dipolemat_domains, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);
//            CDabsorption_domains_Get_Contributions(evecs, evals, dipolemat, parameter_vec, CD, alpha, CIm, GReals, GImags, in, out, p);
            ACD(evecs, evals, Chl_Vec, dipolemat, parameter_vec, ACD_Edge, ACD_Face, CD, alpha, CIm, GReals, GImags, in, out, p);
            get_ACD_RotStrength(evecs, evals, Chl_Vec, dipolemat, RotStrength_vec);

            Linear_Dichroism(evecs, evals, Chl_Vec, dipolemat, parameter_vec, LD, CIm, GReals, GImags, in, out, p);
            get_LD_RotStrength(evecs, evals, Chl_Vec, dipolemat, LDStrength_vec);


            #pragma omp critical
            {
                CD_sum_all(i) = CD(0);
                ACD_Edge_sum_all(i) = ACD_Edge(0);
                ACD_Face_sum_all(i) = ACD_Face(0);
                alpha_sum_all(i) = alpha(0);
                LD_sum_all(i) = LD(0);
                //-----Low Energy (LE) Contributions-----

                CD_sum_LE(i) = CD(1);
                ACD_Edge_sum_LE(i) = ACD_Edge(1);
                ACD_Face_sum_LE(i) = ACD_Face(1);
                alpha_sum_LE(i) = alpha(1);
                LD_sum_LE(i) = LD(1);
                //-----High Energy (HE) Contributions-----

                CD_sum_HE(i) = CD(2);
                ACD_Edge_sum_HE(i) = ACD_Edge(2);
                ACD_Face_sum_HE(i) = ACD_Face(2);
                alpha_sum_HE(i) = alpha(2);
                LD_sum_HE(i) = LD(2);

                RotStrength_sum(i) = RotStrength_vec;
                LDStrength_sum(i) = LDStrength_vec;
            }
            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
        }

        sitevec_random = sitevec;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }

//    FILE * RotStrength_dis_out;
//    RotStrength_dis_out = fopen ("RotStrength.out", "a+");
//
//    fprintf(RotStrength_dis_out,"<R->dis = %lf at Angle %lf \n", dis_RotStrength(RotStrength_sum, -1), t_ang);
//    cout << "dis_RotStrength: " << dis_RotStrength(RotStrength_sum, -1) << endl;
//
//    fclose(RotStrength_dis_out);

    cout << "Loops done! Summing up before writing out..." << endl;

    Vektor<double> CD_sum_total;
    Vektor<double> CD_sum_total_LE;
    Vektor<double> CD_sum_total_HE;

    Vektor<double> ACD_Edge_sum_total;
    Vektor<double> ACD_Edge_sum_total_LE;
    Vektor<double> ACD_Edge_sum_total_HE;

    Vektor<double> ACD_Face_sum_total;
    Vektor<double> ACD_Face_sum_total_LE;
    Vektor<double> ACD_Face_sum_total_HE;

    Vektor<double> alpha_sum_total;
    Vektor<double> alpha_sum_total_LE;
    Vektor<double> alpha_sum_total_HE;

    Vektor<double> LD_sum_total;
    Vektor<double> LD_sum_total_LE;
    Vektor<double> LD_sum_total_HE;

    Vektor< Vektor<double> > RotStrength_sum_total;
    Vektor< Vektor<double> > LDStrength_sum_total;


    CD_sum_total.resize(CD_sum_all(0).size());
    ACD_Edge_sum_total.resize(ACD_Edge_sum_all(0).size());
    ACD_Face_sum_total.resize(ACD_Face_sum_all(0).size());
    alpha_sum_total.resize(alpha_sum_all(0).size());
    LD_sum_total.resize(LD_sum_all(0).size());

    CD_sum_total_LE.resize(CD_sum_LE(0).size());
    ACD_Edge_sum_total_LE.resize(ACD_Edge_sum_LE(0).size());
    ACD_Face_sum_total_LE.resize(ACD_Face_sum_LE(0).size());
    alpha_sum_total_LE.resize(alpha_sum_LE(0).size());
    LD_sum_total_LE.resize(LD_sum_LE(0).size());

    CD_sum_total_HE.resize(CD_sum_HE(0).size());
    ACD_Edge_sum_total_HE.resize(ACD_Edge_sum_HE(0).size());
    ACD_Face_sum_total_HE.resize(ACD_Face_sum_HE(0).size());
    alpha_sum_total_HE.resize(alpha_sum_HE(0).size());
    LD_sum_total_HE.resize(LD_sum_HE(0).size());

    RotStrength_sum_total.resize(RotStrength_sum(0).size());
    LDStrength_sum_total.resize(LDStrength_sum(0).size());

    for(int i = 0; i < (int)CD_sum_all.size(); i++){

        CD_sum_total = CD_sum_total + CD_sum_all(i);
        ACD_Edge_sum_total = ACD_Edge_sum_total + ACD_Edge_sum_all(i);
        ACD_Face_sum_total = ACD_Face_sum_total + ACD_Face_sum_all(i);
        alpha_sum_total = alpha_sum_total + alpha_sum_all(i);
        LD_sum_total = LD_sum_total + LD_sum_all(i);

        CD_sum_total_LE = CD_sum_total_LE + CD_sum_LE(i);
        ACD_Edge_sum_total_LE = ACD_Edge_sum_total_LE + ACD_Edge_sum_LE(i);
        ACD_Face_sum_total_LE = ACD_Face_sum_total_LE + ACD_Face_sum_LE(i);
        alpha_sum_total_LE = alpha_sum_total_LE + alpha_sum_LE(i);
        LD_sum_total_LE = LD_sum_total_LE + LD_sum_LE(i);

        CD_sum_total_HE = CD_sum_total_HE + CD_sum_HE(i);
        ACD_Edge_sum_total_HE = ACD_Edge_sum_total_HE + ACD_Edge_sum_HE(i);
        ACD_Face_sum_total_HE = ACD_Face_sum_total_HE + ACD_Face_sum_HE(i);
        alpha_sum_total_HE = alpha_sum_total_HE + alpha_sum_HE(i);
        LD_sum_total_HE = LD_sum_total_HE + LD_sum_HE(i);

    }

//    cout << "RotStrength_total Size: " << RotStrength_sum_total.size() << endl;
//    cout << "RotStrength Sum(0) : " << RotStrength_sum(0) << endl;
//    cout << "RotStrength Sum(1) : " << RotStrength_sum(1) << endl;
//
//    RotStrength_sum_total = RotStrength_sum(0);
//    RotStrength_sum_total = RotStrength_sum_total + RotStrength_sum(0);

    for(int i = 0; i < (int)RotStrength_sum_total.size(); i++){
    	RotStrength_sum_total(i).resize( (int)RotStrength_sum(0)(i).size() );
    }

    for(int i = 0; i < (int)RotStrength_sum.size(); i++){
    	RotStrength_sum_total += RotStrength_sum(i);
    }

    //-----LD Strengths Summations-----
    for (int i = 0; i < (int) LDStrength_sum_total.size(); i++) {
		LDStrength_sum_total(i).resize((int) LDStrength_sum(0)(i).size());
	}

	for (int i = 0; i < (int) LDStrength_sum.size(); i++) {
		LDStrength_sum_total += LDStrength_sum(i);
	}

//    cout << RotStrength_sum_total << endl;

    alpha_sum_total = alpha_sum_total / parameter_vec(3);
    CD_sum_total = CD_sum_total / parameter_vec(3);
    ACD_Edge_sum_total = ACD_Edge_sum_total / parameter_vec(3);
    ACD_Face_sum_total = ACD_Face_sum_total / parameter_vec(3);
    LD_sum_total = LD_sum_total / parameter_vec(3);

    alpha_sum_total_LE = alpha_sum_total_LE / parameter_vec(3);
    CD_sum_total_LE = CD_sum_total_LE / parameter_vec(3);
    ACD_Edge_sum_total_LE = ACD_Edge_sum_total_LE / parameter_vec(3);
    ACD_Face_sum_total_LE = ACD_Face_sum_total_LE / parameter_vec(3);
    LD_sum_total_LE = LD_sum_total_LE / parameter_vec(3);

    alpha_sum_total_HE = alpha_sum_total_HE / parameter_vec(3);
    CD_sum_total_HE = CD_sum_total_HE / parameter_vec(3);
    ACD_Edge_sum_total_HE = ACD_Edge_sum_total_HE / parameter_vec(3);
    ACD_Face_sum_total_HE = ACD_Face_sum_total_HE / parameter_vec(3);
    LD_sum_total_HE = LD_sum_total_HE / parameter_vec(3);

    for(int i = 0; i < (int)RotStrength_sum_total.size(); i++){
    	RotStrength_sum_total(i) =  RotStrength_sum_total(i) / parameter_vec(3);
    }

    for(int i = 0; i < (int)LDStrength_sum_total.size(); i++){
    	LDStrength_sum_total(i) =  LDStrength_sum_total(i) / parameter_vec(3);
    }


    cout<<"Calculations done, writing out files!"<<endl;

    FILE * ODdisorder;
    FILE * CDdisorder;
    FILE * ACD_Edge_disorder;
    FILE * ACD_Face_disorder;
    FILE * RotStrength_disorder;
    FILE * LDStrength_disorder;
    FILE * LDdisorder;


    string ODFileName = "OD_Dis.out";
    string CDFileName = "CD_Dis.out";
    string ACD_Edge_FileName = "ACD_Edge_Dis.out";
    string ACD_Face_FileName = "ACD_Face_Dis.out";
    string RotStrength_FileName = "RotStrength_Dis.out";
    string LDStrength_FileName = "LDStrength_Dis.out";
    string LDFileName = "LD_Dis.out";

    ODdisorder = fopen ( ODFileName.c_str() , "wb");
    CDdisorder = fopen ( CDFileName.c_str() , "wb");
    ACD_Edge_disorder = fopen ( ACD_Edge_FileName.c_str() , "wb");
    ACD_Face_disorder = fopen ( ACD_Face_FileName.c_str() , "wb");
    RotStrength_disorder = fopen ( RotStrength_FileName.c_str() , "wb");
    LDStrength_disorder = fopen ( LDStrength_FileName.c_str() , "wb");
    LDdisorder = fopen ( LDFileName.c_str() , "wb");

    for(int u = 0; u < (int)CD_sum_total.size(); u++){

		fprintf(ODdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), alpha_sum_total(u), alpha_sum_total_LE(u), alpha_sum_total_HE(u));  //16
		fprintf(CDdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), CD_sum_total(u), CD_sum_total_LE(u), CD_sum_total_HE(u));     //16

		fprintf(ACD_Edge_disorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Edge_sum_total(u), ACD_Edge_sum_total_LE(u), ACD_Edge_sum_total_HE(u));     //16
		fprintf(ACD_Face_disorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Face_sum_total(u), ACD_Face_sum_total_LE(u), ACD_Face_sum_total_HE(u));     //16

		fprintf(LDdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), LD_sum_total(u), LD_sum_total_LE(u), LD_sum_total_HE(u));     //16
	}

    cout << "Printing Rotational Strengths" << endl;
    cout << "RotStrength Sum Total : " << RotStrength_sum_total << endl;


    fprintf(RotStrength_disorder, "Low-Energy Rotational Strengths\n");
    fprintf(RotStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", RotStrength_sum_total(0)(1), RotStrength_sum_total(1)(1), RotStrength_sum_total(2)(1), RotStrength_sum_total(3)(1), RotStrength_sum_total(4)(1), RotStrength_sum_total(5)(1), RotStrength_sum_total(6)(1), RotStrength_sum_total(7)(1), RotStrength_sum_total(8)(1));     //16

    fprintf(RotStrength_disorder, "High-Energy Rotational Strengths\n");
    fprintf(RotStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", RotStrength_sum_total(0)(0), RotStrength_sum_total(1)(0), RotStrength_sum_total(2)(0), RotStrength_sum_total(3)(0), RotStrength_sum_total(4)(0), RotStrength_sum_total(5)(0), RotStrength_sum_total(6)(0), RotStrength_sum_total(7)(0), RotStrength_sum_total(8)(0) );     //16

	cout << "Printing LD Strengths" << endl;
	cout << "LD Sum Total : " << LDStrength_sum_total << endl;

	fprintf(LDStrength_disorder, "Low-Energy Rotational Strengths\n");
	fprintf(LDStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", LDStrength_sum_total(0)(1), LDStrength_sum_total(1)(1), LDStrength_sum_total(2)(1), LDStrength_sum_total(3)(1), LDStrength_sum_total(4)(1), LDStrength_sum_total(5)(1), LDStrength_sum_total(6)(1), LDStrength_sum_total(7)(1), LDStrength_sum_total(8)(1));     //16

	fprintf(LDStrength_disorder, "High-Energy Rotational Strengths\n");
	fprintf(LDStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", LDStrength_sum_total(0)(0), LDStrength_sum_total(1)(0), LDStrength_sum_total(2)(0), LDStrength_sum_total(3)(0), LDStrength_sum_total(4)(0), LDStrength_sum_total(5)(0), LDStrength_sum_total(6)(0), LDStrength_sum_total(7)(0), LDStrength_sum_total(8)(0));     //16



    fclose (ODdisorder);
    fclose (CDdisorder);
    fclose (ACD_Edge_disorder);
    fclose (ACD_Face_disorder);
    fclose (RotStrength_disorder);
    fclose (LDStrength_disorder);
    fclose (LDdisorder);

    cout<<"Writing files done!"<<endl;

}

void ACD_domains(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat, double t_ang, Vektor< ChlAtoms > Chl_Vec){
    //Static Disorder

	cout << "This is ACD_domains() reporting!" << endl;

	Vektor<Vektor<double> > alpha_sum_all((int) parameter_vec(3));
	Vektor<Vektor<double> > CD_sum_all((int) parameter_vec(3));
	Vektor<Vektor<double> > ACD_Edge_sum_all((int) parameter_vec(3));
	Vektor<Vektor<double> > ACD_Face_sum_all((int) parameter_vec(3));
	Vektor<Vektor<double> > LD_sum_all((int) parameter_vec(3));

    Vektor< Vektor< Vektor<double> > > RotStrength_sum ((int)parameter_vec(3));
    Vektor< Vektor< Vektor<double> > > ACD_RotStrength_sum ((int)parameter_vec(3));
    Vektor< Vektor< Vektor<double> > > LDStrength_sum ((int)parameter_vec(3));

    Vektor< Vektor< Vektor<double> > > wM0_sum ((int)parameter_vec(3));

    Vektor< Vektor < Vektor <double> > > Rot_Beitrag_sum((int)parameter_vec(3));
    Vektor< Vektor < Vektor <double> > > ACDRot_Beitrag_sum((int)parameter_vec(3));


	Vektor<Matrix<double> > excmatrix_domains;
	Vektor<Vektor<double> > domains;

	//Ausgliederung von CIm & G(t)
	Vektor<double> CIm;
	correlation(parameter_vec, CIm);

	Vektor<double> GReals;
	Vektor<double> GImags;
	FT_Gt_WZ(GReals, GImags, 8192.0 / 2048.0, 8192.0 * 4.0, parameter_vec);    //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0

	Vektor<Vektor<double> > coupling_vec;
//    build_predefined_excitondomains(chlname_vec, excmatrix_domains, domains, coupling_vec);
	buildexcitondomains(chlname_vec, excmatrix_domains, domains, parameter_vec);

    double ntime = 8192.0 * 4.0;

    printf("Initialization Done! Doing Static-Disorder Loop now.\n");

    int i = 0;
//   #pragma omp parallel for default(none) private(i) shared(dipolemat, ntime,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)
    //	#pragma omp parallel for default(none) private(i) shared(Chl_Vec, dipolemat, ntime,excmatrix_domains, alpha_sum_all, CD_sum_all, ACD_Edge_sum_all, ACD_Face_sum_all, LD_sum_all, parameter_vec,domains,chlname_vec,GReals,GImags,CIm,sitevec)

//    #pragma omp parallel for
    for(i = 0; i < (int)parameter_vec(3); i++){

    	if(i % 50 == 0){
			cout << "Step #: " << i << endl;
    	}

        Vektor< Vektor<double> > CD(3);
        Vektor< Vektor<double> > ACD_Edge(3);
        Vektor< Vektor<double> > ACD_Face(3);
        Vektor< Vektor<double> > alpha(3);
        Vektor< Vektor<double> > LD(3);

        Matrix<double> evecs;            //matrix&vector for Eigenvectors & Eigenvalues
        Vektor<double> evals;

        Vektor< Vektor<double> > RotStrength_vec( excmatrix_domains.size() );
        Vektor< Vektor<double> > ACD_RotStrength_vec( excmatrix_domains.size() );
        Vektor< Vektor<double> > LDStrength_vec( excmatrix_domains.size() );
        Vektor< Vektor<double> > wM0_vec( excmatrix_domains.size() );

        Vektor < Vektor <double> > Rot_Beitrag;
        Vektor < Vektor <double> > ACDRot_Beitrag;

        Vektor<double> sitevec_random;
        sitevec_random = sitevec;

        fftw_complex *in, *out;
        fftw_plan p;

		in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);     //these vectors are empty, contrain 0's!
		out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * ntime);

		RandomizeSiteEnergies(sitevec_random, parameter_vec(2) / 2.35482);
//        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
		addsiteenergies(excmatrix_domains, domains, sitevec_random);

        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;

        for(int j = 0; j < (int)excmatrix_domains.size(); j++){  //(int)excmatrix_domains.size()

			Vektor<Vektor<Vektor<double> > > dipolemat_domains;
			dipolemat_domains = dipolemat;

			//cout<<"Building dipolemat_domains :  "<<endl;

			dipolemat_domains(0).resize(domains(j).size());
			dipolemat_domains(1).resize(domains(j).size());

			for (int u = 0; u < (int) dipolemat_domains(0).size(); u++) {

				dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u) - 1);
				dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u) - 1);
			}

			Vektor< ChlAtoms > Chl_Vec_domains;
			Chl_Vec_domains.resize(domains(j).size());

			for(int u = 0; u < (int) domains(j).size(); u++){
				Chl_Vec_domains(u) = Chl_Vec(domains(j)(u) - 1);
			}

//			for(int u = 0; u < (int)Chl_Vec_domains.size(); u++){
//				cout << Chl_Vec_domains(u).Dipole_Mom << "  ";
//			}
//			cout << endl;
////			cout << endl << Chl_Vec_domains.size() << endl;
////			cout << Chl_Vec.size() << endl;
//			cout << dipolemat_domains << endl;

            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));

            ACD(evecs, evals, Chl_Vec_domains, dipolemat_domains, parameter_vec, ACD_Edge, ACD_Face, CD, alpha, CIm, GReals, GImags, in, out, p);
            Linear_Dichroism(evecs, evals, Chl_Vec_domains, dipolemat_domains, parameter_vec, LD, CIm, GReals, GImags, in, out, p);

//            get_RotStrengths(evecs, evals, Chl_Vec_domains, dipolemat_domains, LDStrength_vec(j), RotStrength_vec(j), ACD_RotStrength_vec(j), CIm, parameter_vec, wM0_vec(j), Rot_Beitrag, ACDRot_Beitrag);

//            //-----Output of Rotstrengths-----
//
//			ofstream RotStrength_dis_out;
//			string file_name = "RotStrength_" + std::to_string(j);
//			RotStrength_dis_out.open(file_name);
//
//			RotStrength_dis_out << evals << endl << RotStrength_vec(j) << endl << endl;
//			RotStrength_dis_out << evals << endl << ACD_RotStrength_vec(j) << endl << endl;
//
//			RotStrength_dis_out.close();



//            get_LD_RotStrength(evecs, evals, Chl_Vec, dipolemat, LDStrength_vec);
//            get_ACD_RotStrength(evecs, evals, Chl_Vec, dipolemat, RotStrength_vec);

//            #pragma omp critical
            {
                CD_sum_all(i) = CD(0);
                ACD_Edge_sum_all(i) = ACD_Edge(0);
                ACD_Face_sum_all(i) = ACD_Face(0);
                alpha_sum_all(i) = alpha(0);
                LD_sum_all(i) = LD(0);

                RotStrength_sum(i) = RotStrength_vec;
                LDStrength_sum(i) = LDStrength_vec;
                ACD_RotStrength_sum(i) = ACD_RotStrength_vec;
                wM0_sum(i) = wM0_vec;

                Rot_Beitrag_sum(i) = Rot_Beitrag;
                ACDRot_Beitrag_sum(i) = ACDRot_Beitrag;
            }
            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
        }

        sitevec_random = sitevec;
        fftw_destroy_plan(p);
        fftw_free(in);
        fftw_free(out);
    }


    cout << "Loops done! Summing up before writing out..." << endl;


    Vektor<double> CD_sum_total;
    Vektor<double> ACD_Edge_sum_total;
    Vektor<double> ACD_Face_sum_total;
    Vektor<double> alpha_sum_total;
    Vektor<double> LD_sum_total;

    CD_sum_total.resize(CD_sum_all(0).size());
    ACD_Edge_sum_total.resize(ACD_Edge_sum_all(0).size());
    ACD_Face_sum_total.resize(ACD_Face_sum_all(0).size());
    alpha_sum_total.resize(alpha_sum_all(0).size());
    LD_sum_total.resize(LD_sum_all(0).size());


    //-----Sum up the Rotational / ACD / LD Strengths-----
    Vektor<Vektor<double> > RotStrength_sum_total;
    Vektor<Vektor<double> > ACD_RotStrength_sum_total;
	Vektor<Vektor<double> > LDStrength_sum_total;
	Vektor<Vektor<double> > wM0_sum_total;

	Vektor < Vektor <double> >  Rot_Beitrag_sum_total;
	Vektor < Vektor <double> >  ACDRot_Beitrag_sum_total;


	RotStrength_sum_total.resize(RotStrength_sum(0).size());
	ACD_RotStrength_sum_total.resize(ACD_RotStrength_sum(0).size());
	LDStrength_sum_total.resize(LDStrength_sum(0).size());
	wM0_sum_total.resize(wM0_sum(0).size());

	Rot_Beitrag_sum_total.resize(Rot_Beitrag_sum(0).size());
	ACDRot_Beitrag_sum_total.resize(ACDRot_Beitrag_sum(0).size());

//-----This is from the LE & HE Contributions-----
//    Vektor<double> CD_sum_total_LE;
//    Vektor<double> CD_sum_total_HE;
//
//    Vektor<double> ACD_Edge_sum_total_LE;
//    Vektor<double> ACD_Edge_sum_total_HE;
//
//    Vektor<double> ACD_Face_sum_total_LE;
//    Vektor<double> ACD_Face_sum_total_HE;
//
//    Vektor<double> alpha_sum_total_LE;
//    Vektor<double> alpha_sum_total_HE;
//
//    Vektor<double> LD_sum_total_LE;
//    Vektor<double> LD_sum_total_HE;
//

//
//    CD_sum_total_LE.resize(CD_sum_LE(0).size());
//    ACD_Edge_sum_total_LE.resize(ACD_Edge_sum_LE(0).size());
//    ACD_Face_sum_total_LE.resize(ACD_Face_sum_LE(0).size());
//    alpha_sum_total_LE.resize(alpha_sum_LE(0).size());
//    LD_sum_total_LE.resize(LD_sum_LE(0).size());
//
//    CD_sum_total_HE.resize(CD_sum_HE(0).size());
//    ACD_Edge_sum_total_HE.resize(ACD_Edge_sum_HE(0).size());
//    ACD_Face_sum_total_HE.resize(ACD_Face_sum_HE(0).size());
//    alpha_sum_total_HE.resize(alpha_sum_HE(0).size());
//    LD_sum_total_HE.resize(LD_sum_HE(0).size());
//
//    RotStrength_sum_total.resize(RotStrength_sum(0).size());
//    LDStrength_sum_total.resize(LDStrength_sum(0).size());
//


//    cout << ACD_Edge_sum_all.size() << endl;
//    cout << CD_sum_all.size() << endl;
//    cout << ACD_Edge_sum_all.size() << endl;
//    cout << ACD_Face_sum_all.size() << endl;
//    cout << alpha_sum_all.size() << endl;
//    cout << LD_sum_all.size() << endl;
//
//    exit(1);

    for(int i = 0; i < (int)CD_sum_all.size(); i++){

        CD_sum_total = CD_sum_total + CD_sum_all(i);
        ACD_Edge_sum_total = ACD_Edge_sum_total + ACD_Edge_sum_all(i);
        ACD_Face_sum_total = ACD_Face_sum_total + ACD_Face_sum_all(i);
        alpha_sum_total = alpha_sum_total + alpha_sum_all(i);
        LD_sum_total = LD_sum_total + LD_sum_all(i);

//        CD_sum_total_LE = CD_sum_total_LE + CD_sum_LE(i);
//        ACD_Edge_sum_total_LE = ACD_Edge_sum_total_LE + ACD_Edge_sum_LE(i);
//        ACD_Face_sum_total_LE = ACD_Face_sum_total_LE + ACD_Face_sum_LE(i);
//        alpha_sum_total_LE = alpha_sum_total_LE + alpha_sum_LE(i);
//        LD_sum_total_LE = LD_sum_total_LE + LD_sum_LE(i);
//
//        CD_sum_total_HE = CD_sum_total_HE + CD_sum_HE(i);
//        ACD_Edge_sum_total_HE = ACD_Edge_sum_total_HE + ACD_Edge_sum_HE(i);
//        ACD_Face_sum_total_HE = ACD_Face_sum_total_HE + ACD_Face_sum_HE(i);
//        alpha_sum_total_HE = alpha_sum_total_HE + alpha_sum_HE(i);
//        LD_sum_total_HE = LD_sum_total_HE + LD_sum_HE(i);

    }

    //-----Sum up the Rotational / ACD / LD Strengths-----
	for (int i = 0; i < (int) RotStrength_sum_total.size(); i++) {
		RotStrength_sum_total(i).resize((int) RotStrength_sum(0)(i).size());
		ACD_RotStrength_sum_total(i).resize((int) ACD_RotStrength_sum(0)(i).size());
		LDStrength_sum_total(i).resize((int) LDStrength_sum(0)(i).size());
		wM0_sum_total(i).resize((int) wM0_sum(0)(i).size());
	}

	for (int i = 0; i < (int) RotStrength_sum.size(); i++) {
		RotStrength_sum_total += RotStrength_sum(i);
		ACD_RotStrength_sum_total += ACD_RotStrength_sum(i);
		LDStrength_sum_total += LDStrength_sum(i);
		wM0_sum_total += wM0_sum(i);
	}

	//-----Normalize the Strengths-----
	for(int i = 0; i < (int)  RotStrength_sum_total.size(); i++) {
		RotStrength_sum_total(i) = RotStrength_sum_total(i) / parameter_vec(3);
		ACD_RotStrength_sum_total(i) = ACD_RotStrength_sum_total(i) / parameter_vec(3);
		LDStrength_sum_total(i) = LDStrength_sum_total(i) / parameter_vec(3);
		wM0_sum_total(i) = wM0_sum_total(i) / parameter_vec(3);
	}

	//-----Output the Inh. Strengths-----
	cout << RotStrength_sum_total << endl;
	cout << ACD_RotStrength_sum_total << endl;
	cout << LDStrength_sum_total << endl;
	cout << wM0_sum_total << endl;

	//-----Rotational Strength Beitrag-----
	for(int i = 0; i < (int) Rot_Beitrag_sum_total.size(); i++){
		Rot_Beitrag_sum_total(i).resize((int) Rot_Beitrag_sum(0)(i).size());
	}

	for(int i = 0; i < (int) Rot_Beitrag_sum.size(); i++){
		Rot_Beitrag_sum_total += Rot_Beitrag_sum(i);
	}

	for(int i = 0; i < (int) Rot_Beitrag_sum_total.size(); i++){
		Rot_Beitrag_sum_total(i) = Rot_Beitrag_sum_total(i) / parameter_vec(3);
	}

	ofstream Out_RotBeitrag;
	Out_RotBeitrag.open("RM_Strength_Beitrag.out");

	for(int i = 0; i < (int)Rot_Beitrag_sum_total.size(); i++){

		for(int j = 0; j < (int) Rot_Beitrag_sum_total(i).size(); j++){
			Out_RotBeitrag << Rot_Beitrag_sum_total(i)(j) << "  ";
		}

		Out_RotBeitrag << endl;
//		Out_RotBeitrag << Rot_Beitrag_sum_total(i) << endl;
	}

	Out_RotBeitrag.close();

	//-----ACD Rotational Strength Beitrag-----
	for(int i = 0; i < (int) ACDRot_Beitrag_sum_total.size(); i++){
		ACDRot_Beitrag_sum_total(i).resize((int) ACDRot_Beitrag_sum(0)(i).size());
	}

	for(int i = 0; i < (int) ACDRot_Beitrag_sum.size(); i++){
		ACDRot_Beitrag_sum_total += ACDRot_Beitrag_sum(i);
	}

	for(int i = 0; i < (int) ACDRot_Beitrag_sum_total.size(); i++){
		ACDRot_Beitrag_sum_total(i) = ACDRot_Beitrag_sum_total(i) / parameter_vec(3);
	}

	ofstream Out_ACDRotBeitrag;
	Out_ACDRotBeitrag.open("RM_ACDStrength_Beitrag.out");

	for(int i = 0; i < (int)ACDRot_Beitrag_sum_total.size(); i++){

		for(int j = 0; j < (int) ACDRot_Beitrag_sum_total(i).size(); j++){
			Out_ACDRotBeitrag << ACDRot_Beitrag_sum_total(i)(j) << "  ";
		}

		Out_ACDRotBeitrag << endl;
//		Out_RotBeitrag << Rot_Beitrag_sum_total(i) << endl;
	}

	Out_ACDRotBeitrag.close();


////    cout << "RotStrength_total Size: " << RotStrength_sum_total.size() << endl;
////    cout << "RotStrength Sum(0) : " << RotStrength_sum(0) << endl;
////    cout << "RotStrength Sum(1) : " << RotStrength_sum(1) << endl;
////
////    RotStrength_sum_total = RotStrength_sum(0);
////    RotStrength_sum_total = RotStrength_sum_total + RotStrength_sum(0);
//
//    for(int i = 0; i < (int)RotStrength_sum_total.size(); i++){
//    	RotStrength_sum_total(i).resize( (int)RotStrength_sum(0)(i).size() );
//    }
//
//    for(int i = 0; i < (int)RotStrength_sum.size(); i++){
//    	RotStrength_sum_total += RotStrength_sum(i);
//    }
//
//    //-----LD Strengths Summations-----
//    for (int i = 0; i < (int) LDStrength_sum_total.size(); i++) {
//		LDStrength_sum_total(i).resize((int) LDStrength_sum(0)(i).size());
//	}
//
//	for (int i = 0; i < (int) LDStrength_sum.size(); i++) {
//		LDStrength_sum_total += LDStrength_sum(i);
//	}
//
////    cout << RotStrength_sum_total << endl;
//
    alpha_sum_total = alpha_sum_total / parameter_vec(3);
    CD_sum_total = CD_sum_total / parameter_vec(3);
    ACD_Edge_sum_total = ACD_Edge_sum_total / parameter_vec(3);
    ACD_Face_sum_total = ACD_Face_sum_total / parameter_vec(3);
    LD_sum_total = LD_sum_total / parameter_vec(3);

//    alpha_sum_total_LE = alpha_sum_total_LE / parameter_vec(3);
//    CD_sum_total_LE = CD_sum_total_LE / parameter_vec(3);
//    ACD_Edge_sum_total_LE = ACD_Edge_sum_total_LE / parameter_vec(3);
//    ACD_Face_sum_total_LE = ACD_Face_sum_total_LE / parameter_vec(3);
//    LD_sum_total_LE = LD_sum_total_LE / parameter_vec(3);
//
//    alpha_sum_total_HE = alpha_sum_total_HE / parameter_vec(3);
//    CD_sum_total_HE = CD_sum_total_HE / parameter_vec(3);
//    ACD_Edge_sum_total_HE = ACD_Edge_sum_total_HE / parameter_vec(3);
//    ACD_Face_sum_total_HE = ACD_Face_sum_total_HE / parameter_vec(3);
//    LD_sum_total_HE = LD_sum_total_HE / parameter_vec(3);
//
//
//
    cout<<"Calculations done, writing out files!"<<endl;

    FILE * ODdisorder;
    FILE * CDdisorder;
    FILE * ACD_Edge_disorder;
    FILE * ACD_Face_disorder;
    FILE * LDdisorder;

//    FILE * RotStrength_disorder;
//    FILE * LDStrength_disorder;
//
//
    string ODFileName = "OD_Dis.out";
    string CDFileName = "CD_Dis.out";
    string ACD_Edge_FileName = "ACD_Edge_Dis.out";
    string ACD_Face_FileName = "ACD_Face_Dis.out";
    string LDFileName = "LD_Dis.out";

//    string RotStrength_FileName = "RotStrength_Dis.out";
//    string LDStrength_FileName = "LDStrength_Dis.out";


    ODdisorder = fopen ( ODFileName.c_str() , "wb");
    CDdisorder = fopen ( CDFileName.c_str() , "wb");
    ACD_Edge_disorder = fopen ( ACD_Edge_FileName.c_str() , "wb");
    ACD_Face_disorder = fopen ( ACD_Face_FileName.c_str() , "wb");
    LDdisorder = fopen ( LDFileName.c_str() , "wb");

//    RotStrength_disorder = fopen ( RotStrength_FileName.c_str() , "wb");
//    LDStrength_disorder = fopen ( LDStrength_FileName.c_str() , "wb");


    for (int u = 0; u < (int) CD_sum_total.size(); u++) {

		fprintf(ODdisorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), alpha_sum_total(u) );  //16
		fprintf(CDdisorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), CD_sum_total(u) );     //16

		fprintf(ACD_Edge_disorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Edge_sum_total(u) );     //16
		fprintf(ACD_Face_disorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Face_sum_total(u) );     //16

		fprintf(LDdisorder, " %f  %f \n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), LD_sum_total(u) );     //16
	}


//    for(int u = 0; u < (int)CD_sum_total.size(); u++){
//
//		fprintf(ODdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), alpha_sum_total(u), alpha_sum_total_LE(u), alpha_sum_total_HE(u));  //16
//		fprintf(CDdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), CD_sum_total(u), CD_sum_total_LE(u), CD_sum_total_HE(u));     //16
//
//		fprintf(ACD_Edge_disorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Edge_sum_total(u), ACD_Edge_sum_total_LE(u), ACD_Edge_sum_total_HE(u));     //16
//		fprintf(ACD_Face_disorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), ACD_Face_sum_total(u), ACD_Face_sum_total_LE(u), ACD_Face_sum_total_HE(u));     //16
//
//		fprintf(LDdisorder, " %f  %f  %f  %f\n", 1e9 * (0.01 * 1.0 / (1.0 / (300 * 1e-7) - u * (2 * 3.1415926) / (4.0))), LD_sum_total(u), LD_sum_total_LE(u), LD_sum_total_HE(u));     //16
//	}

//    cout << "Printing Rotational Strengths" << endl;
//    cout << "RotStrength Sum Total : " << RotStrength_sum_total << endl;
//
//
//    fprintf(RotStrength_disorder, "Low-Energy Rotational Strengths\n");
//    fprintf(RotStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", RotStrength_sum_total(0)(1), RotStrength_sum_total(1)(1), RotStrength_sum_total(2)(1), RotStrength_sum_total(3)(1), RotStrength_sum_total(4)(1), RotStrength_sum_total(5)(1), RotStrength_sum_total(6)(1), RotStrength_sum_total(7)(1), RotStrength_sum_total(8)(1));     //16
//
//    fprintf(RotStrength_disorder, "High-Energy Rotational Strengths\n");
//    fprintf(RotStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", RotStrength_sum_total(0)(0), RotStrength_sum_total(1)(0), RotStrength_sum_total(2)(0), RotStrength_sum_total(3)(0), RotStrength_sum_total(4)(0), RotStrength_sum_total(5)(0), RotStrength_sum_total(6)(0), RotStrength_sum_total(7)(0), RotStrength_sum_total(8)(0) );     //16
//
//	cout << "Printing LD Strengths" << endl;
//	cout << "LD Sum Total : " << LDStrength_sum_total << endl;
//
//	fprintf(LDStrength_disorder, "Low-Energy Rotational Strengths\n");
//	fprintf(LDStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", LDStrength_sum_total(0)(1), LDStrength_sum_total(1)(1), LDStrength_sum_total(2)(1), LDStrength_sum_total(3)(1), LDStrength_sum_total(4)(1), LDStrength_sum_total(5)(1), LDStrength_sum_total(6)(1), LDStrength_sum_total(7)(1), LDStrength_sum_total(8)(1));     //16
//
//	fprintf(LDStrength_disorder, "High-Energy Rotational Strengths\n");
//	fprintf(LDStrength_disorder, " %f  %f  %f  %f  %f  %f  %f  %f  %f\n", LDStrength_sum_total(0)(0), LDStrength_sum_total(1)(0), LDStrength_sum_total(2)(0), LDStrength_sum_total(3)(0), LDStrength_sum_total(4)(0), LDStrength_sum_total(5)(0), LDStrength_sum_total(6)(0), LDStrength_sum_total(7)(0), LDStrength_sum_total(8)(0));     //16
//
//
//
    fclose (ODdisorder);
    fclose (CDdisorder);
    fclose (ACD_Edge_disorder);
    fclose (ACD_Face_disorder);
    fclose (LDdisorder);
//    fclose (RotStrength_disorder);
//    fclose (LDStrength_disorder);
//
    cout<<"Writing files done!"<<endl;

}

//double StaticDisorder_domains_stoerungstheorie(Vektor<double> sitevec, Vektor<double> parameter_vec, Vektor<string> chlname_vec,  Vektor< Vektor< Vektor<double> > > dipolemat){
//
//    Vektor< Vektor<double> > alpha_sum ((int)parameter_vec(3));
//    Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
//
//    Vektor< Matrix<double> > excmatrix_domains;
//    Vektor< Vektor<double> > domains;
//
//    //Ausgliederung von CIm & G(t)
//    Vektor<double> CIm;
//    correlation(parameter_vec, CIm);
//
//    Vektor<double> GReals;
//    Vektor<double> GImags;
//    FT_Gt_WZ(GReals,GImags,8192.0/2048.0 ,8192.0 * 4.0,parameter_vec);    //original: 8192.0 / 16.0 ,8192.0*512 //fast version 8192.0/256.0 ,8192.0 * 32.0  //Normal Version 22.06.2015 8192.0/512.0 ,8192.0 * 16.0
//
//    Vektor< Vektor<double> > coupling_vec;    //holds all couplings from excitonic_couplings.dat
//    build_predefined_excitondomains(chlname_vec, excmatrix_domains, domains, coupling_vec);
//
//    Matrix<double> coupling_mat;
//    MakeCouplingMatrix(coupling_vec, coupling_mat, chlname_vec.size() + 1);
//
////    cout << "GetCouplingsFunction : " << getcouplings(coupling_vec,50,51) << endl;
////    cout << "Coupling_Mat : " << coupling_mat(50,51) <<endl;
////    exit(1);
//
//    //cout<<"Domainentry :  "<<domains(0)(0)<<endl;
//    //cout<<endl<<endl<<"Dipolemat :  "<<dipolemat<<endl;
//    //cout<<"Excmatrix_domains size :  "<<excmatrix_domains.size()<<endl;
//
//    int i = 0;
//    #pragma omp parallel for default(none) private(i) shared(dipolemat,excmatrix_domains,CD_sum,alpha_sum,parameter_vec,domains,chlname_vec,coupling_mat,GReals,GImags,CIm,sitevec)
//    for(i = 0; i < (int)parameter_vec(3); i++){ //parameter_vec(3) holds #nrand
//
////        if(i % 25 == 0){
////
////            cout << "We are currently at iterationstep : " << i <<endl;
////        }
//
//        Vektor<double> CD;
//        Vektor<double> alpha;
//        Matrix<double> evecs;            //matrix&vector for Eigenvectors & Eigenvalues
//        Vektor<double> evals;
//
//        Matrix<double> evecshigher;            //matrix&vector for Eigenvectors & Eigenvalues
//        Vektor<double> evalshigher;
//
//        Vektor<double> sitevec_random;
//        sitevec_random = sitevec;
//
//        RandomizeSiteEnergies(sitevec_random,parameter_vec(2)/2.35482);  //130.0/2.35482 set temporarily //parameter_vec(2)/2.35482 contains FWHM = 130cm^-1//The boost Normaldistribution has the inputform (MEAN-Value, Sigma)
//        addsiteenergies(excmatrix_domains,domains,sitevec_random);
//
//        int stdomain = 4;
//
//        //cout<<"Excmatrix_domains : "<<endl<<excmatrix_domains<<endl;
//
//        for(int j = 0; j < (int)excmatrix_domains.size() - 1; j++){  //(int)excmatrix_domains.size()
//
//            Vektor< Vektor< Vektor<double> > > dipolemat_domains;
//            dipolemat_domains = dipolemat;
//
//
//            int domainsizeQy = domains(j).size();
//            int domainsizeHigher = domains(stdomain).size();
//
//            dipolemat_domains(0).resize(domainsizeQy + domainsizeHigher);
//            dipolemat_domains(1).resize(domainsizeQy + domainsizeHigher);
//
//            for(int u = 0; u < domainsizeQy; u++){
//
//                dipolemat_domains(0)(u) = dipolemat(0)(domains(j)(u)-1);
//                dipolemat_domains(1)(u) = dipolemat(1)(domains(j)(u)-1);
//            }
//
//            for(int u = 0; u < domainsizeHigher; u++){
//
//                dipolemat_domains(0)(u + domainsizeQy) = dipolemat(0)(domains(stdomain)(u)-1);
//                dipolemat_domains(1)(u + domainsizeQy) = dipolemat(1)(domains(stdomain)(u)-1);
//            }
//
//            diagonalizeexcmatrix_domains(evecs, evals, chlname_vec, parameter_vec, excmatrix_domains(j));
//            diagonalizeexcmatrix_domains(evecshigher, evalshigher, chlname_vec, parameter_vec, excmatrix_domains(stdomain));
//
//
////            cout << "Dipolemat that was initiated : " << endl << endl << endl;
////            cout << dipolemat_domains << endl;
////
////            cout<<"Exciton Matrix : "<<endl<<endl;
////            cout << excmatrix_domains(j) <<endl<<endl;
////
////            cout<<"evals: "<<evals<<endl<<endl<<endl;
////            cout<<"evecs: "<<evecs<<endl<<endl<<endl;
////
////            cout<<"evals-higher: "<<evalshigher<<endl<<endl<<endl;
////            cout<<"evecs-higher: "<<evecshigher<<endl<<endl<<endl;
//
//            stoerungstheorie(evecs,evals,evalshigher,evecshigher,coupling_mat,domains,j,stdomain);
//
////            cout << "before Absorption" <<endl;
////
////            cout << "Evals new: "<<endl<<endl<<evals<<endl<<endl;
////
////            cout << "Evecs new: "<<endl<<endl;
////            cout<<evecs<<endl<<endl;
////
//////            for(int i = 0; i < evecs.size2(); i++){
//////                cout << evecs(11,i)<<endl;
//////            }
////
////            cout << "dipolemat domains : " <<  dipolemat_domains(0) <<endl;
////            cout << "dipolemat domains centers : "<< dipolemat_domains(1) <<endl;
////
////            exit(1);
//
//            CDabsorption_domains(evecs, evals, dipolemat_domains, parameter_vec, CD, alpha, CIm, GReals, GImags);
//
//
//            #pragma omp critical
//            {
//                CD_sum(i) = CD;
//                alpha_sum(i) = alpha;
//            }
//
//            //Vektor< Vektor<double> > CD_sum ((int)parameter_vec(3));
//        }
//
//        sitevec_random = sitevec;
//    }
//
//    Vektor<double> CD_sum_total;
//    Vektor<double> alpha_sum_total;
//
//    CD_sum_total.resize(CD_sum(0).size());
//    alpha_sum_total.resize(alpha_sum(0).size());
//
//    for(int i = 0; i < (int)CD_sum.size(); i++){
//        CD_sum_total = CD_sum_total + CD_sum(i);
//        alpha_sum_total = alpha_sum_total + alpha_sum(i);
//    }
//
//    alpha_sum_total = alpha_sum_total / parameter_vec(3);
//    CD_sum_total = CD_sum_total / parameter_vec(3);
//
//    cout<<"Calculations done, writing out files!"<<endl;
//
//    FILE * ODdisorder;
//    FILE * CDdisorder;
//
//    ODdisorder = fopen ("ODDis.out", "wb");
//    CDdisorder = fopen ("CDDis.out", "wb");
//
//    for(int u = 0; u < (int)CD_sum_total.size(); u++){
//
//        fprintf(ODdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),alpha_sum_total(u));  //16
//        fprintf(CDdisorder," %f  %f\n", 1e9 * (0.01 * 1.0/(1.0/(300*1e-7) - u*(2*3.1415926)/(4.0))),CD_sum_total(u));     //16
//    }
//
//    fclose (ODdisorder);
//    fclose (CDdisorder);
//
//    cout<<"Writing files done!"<<endl;
//
//    return 0;
//}

#endif // RANDOMGAUSS_HPP_INCLUDED
