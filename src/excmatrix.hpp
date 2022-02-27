#ifndef EXCMATRIX_HPP_INCLUDED
#define EXCMATRIX_HPP_INCLUDED

#include <iostream>
#include "Utilities.hpp"
#include <fstream>
#include <sstream>
#include "nr3.h"
#include "eigen_sym.h"
#include "eigen_unsym.h"

using namespace std;

void getSiteEnergies(Vektor<double>& sitevec, Vektor<string> chlname_vec){
    //at first get the site energies out of the .dat file and save it to a site-energy-vector sitevec
    ifstream Site_Energies ("src/Site_Energies.dat");

    double E0a,E0b,E0aBy,E0bBy,E0aQx,E0bQx,E0aBx,E0bBx,LUT,LUT01,LUT02,LUT03,NEO,NEO01,NEO02,NEO03,VIO,VIO01,VIO02,VIO03,E0aN,E0bN;                     //may needs to be expaned if carotenoids || Qx,Bx etc. transitions are present
    string str;


    if (Site_Energies.is_open()){           //obtain E0a & E0b from Site_Energies.dat
        if(Site_Energies.good()){           //this is from the first line
            getline (Site_Energies,str);
            stringstream sstest(str);
            sstest >> E0a;
            sstest >> E0b;
            sstest >> E0aBy;
            sstest >> E0bBy;
            sstest >> E0aQx;
            sstest >> E0bQx;
            sstest >> E0aBx;
            sstest >> E0bBx;
            sstest >> LUT;
            sstest >> LUT01;
            sstest >> LUT02;
            sstest >> LUT03;
            sstest >> VIO;
            sstest >> VIO01;
            sstest >> VIO02;
            sstest >> VIO03;
            sstest >> NEO;
            sstest >> NEO01;
            sstest >> NEO02;
            sstest >> NEO03;
            sstest >> E0aN;
            sstest >> E0bN;
        }
    }
    cout<<"E0a: "<<E0a<<"  E0b: "<<E0b<<" E0aBy: "<<E0aBy<<" E0bBy: "<<E0bBy<<" E0aQx: "<<E0aQx<<" E0bQx: "<<E0bQx<<" E0aBx: "<<E0aBx<<" E0bBx: "<<E0bBx<<" LUT: "<<LUT<<" LUT01: "<<LUT01<<" LUT02: "<<LUT02<<" LUT03: "<<LUT03<<" VIO: "<<VIO
    <<" VIO01: "<<VIO01<<" VIO02: "<<VIO02<<" VIO03: "<<VIO03<<" NEO: "<<NEO<<" NEO01: "<<NEO01<<" NEO02: "<<NEO02<<" NEO03: "<<NEO03<<" E0aN: "<<E0aN<<" E0bN: "<<E0bN<<endl;     //errorchecking

    int k=0;
    if (Site_Energies.is_open()){           //obtain siteenergies from Site_Energies.dat
        while(Site_Energies.good()){           //this is from all other lines
            getline (Site_Energies,str);
            stringstream sstest(str);
            sitevec.resize(sitevec.size()+1);
            sstest >> sitevec(k);
            k++;
        }
    }

    //resize the siteenergyvector so that there are no 0's at the end ---- this is an arbitrary and temporary choice!
    while(sitevec(sitevec.size()-1)==0){
        sitevec.resize(sitevec.size()-1);
    }
    //cout<<sitevec<<endl;        //errorchecking

    //assigns E0a & E0b to the site energies
    //cout<<chlname_vec<<endl;

    int n = sitevec.size();
	for (int i = 0; i < n; i++) {
		if (chlname_vec(i) == "CLA") {
			sitevec(i) = sitevec(i) + 1.0 / (E0a * 1e-7);
		}
		if (chlname_vec(i) == "CHL") {
			sitevec(i) = sitevec(i) + 1.0 / (E0b * 1e-7);
		}
		if (chlname_vec(i) == "CLABy") {
			sitevec(i) = sitevec(i) + 1.0 / (E0aBy * 1e-7);
		}
		if (chlname_vec(i) == "CHLBy") {
			sitevec(i) = sitevec(i) + 1.0 / (E0bBy * 1e-7);
		}
		if (chlname_vec(i) == "CLAQx") {
			sitevec(i) = sitevec(i) + 1.0 / (E0aQx * 1e-7);
		}
		if (chlname_vec(i) == "CHLQx") {
			sitevec(i) = sitevec(i) + 1.0 / (E0bQx * 1e-7);
		}
		if (chlname_vec(i) == "CLABx") {
			sitevec(i) = sitevec(i) + 1.0 / (E0aBx * 1e-7);
		}
		if (chlname_vec(i) == "CHLBx") {
			sitevec(i) = sitevec(i) + 1.0 / (E0bBx * 1e-7);
		}
		if (chlname_vec(i) == "LUT") {
			sitevec(i) = sitevec(i) + 1.0 / (LUT * 1e-7);
		}
		if (chlname_vec(i) == "LUT01") {
			sitevec(i) = sitevec(i) + 1.0 / (LUT01 * 1e-7);
		}
		if (chlname_vec(i) == "LUT02") {
			sitevec(i) = sitevec(i) + 1.0 / (LUT02 * 1e-7);
		}
		if (chlname_vec(i) == "LUT03") {
			sitevec(i) = sitevec(i) + 1.0 / (LUT03 * 1e-7);
		}
		if (chlname_vec(i) == "XAT") {
			sitevec(i) = sitevec(i) + 1.0 / (VIO * 1e-7);
		}
		if (chlname_vec(i) == "XAT01") {
			sitevec(i) = sitevec(i) + 1.0 / (VIO01 * 1e-7);
		}
		if (chlname_vec(i) == "XAT02") {
			sitevec(i) = sitevec(i) + 1.0 / (VIO02 * 1e-7);
		}
		if (chlname_vec(i) == "XAT03") {
			sitevec(i) = sitevec(i) + 1.0 / (VIO03 * 1e-7);
		}
		if (chlname_vec(i) == "NEX") {
			sitevec(i) = sitevec(i) + 1.0 / (NEO * 1e-7);
		}
		if (chlname_vec(i) == "NEX01") {
			sitevec(i) = sitevec(i) + 1.0 / (NEO01 * 1e-7);
		}
		if (chlname_vec(i) == "NEX02") {
			sitevec(i) = sitevec(i) + 1.0 / (NEO02 * 1e-7);
		}
		if (chlname_vec(i) == "NEX03") {
			sitevec(i) = sitevec(i) + 1.0 / (NEO03 * 1e-7);
		}
		if (chlname_vec(i) == "CLAN") {
			sitevec(i) = sitevec(i) + 1.0 / (E0aN * 1e-7);
		}
		if (chlname_vec(i) == "CHLN") {
			sitevec(i) = sitevec(i) + 1.0 / (E0bN * 1e-7);
		}

		//Expand for other tranitions!!
	}
	//cout<<"SiteVec is: "<<sitevec<<endl;  //errorchecking
}

void diagonalizeexcmatrix_domains(Matrix<double>& evec, Vektor<double>& eval, Vektor<string>& chlname_vec, Vektor<double>& parameter_vec, Matrix<double>& couplingmatrix){

    string str;

    int n = couplingmatrix.size1();

    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Here the Eigenvalueproblem is solved through the NR eigensolver "nr3.h" && eigen_sym
    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //NR MatrixClass .... this is just a temporary object
    MatDoub tmp_excmat(n,n);

    //Copy the excitonmatrix to the NR-Matrixclass

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tmp_excmat[i][j] = couplingmatrix(i, j);
		}
	}

    //NR eigenvalueproblem solver - this creats a class
    Symmeig sol(tmp_excmat);
    MatDoub tmp_evec(n,n);
    VecDoub tmp_eval(n);
	tmp_eval = sol.d;     //obtain solutions from the solverclass
    tmp_evec = sol.z;

    eval.resize(n);         //resize original working objects
    evec.resize(n,n);

	//Copy Eigenvalues
	for (int i = 0; i < n; i++) {
		eval(i) = tmp_eval[i];
	}

	//Copy Eigencoefficients
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			evec(i, j) = tmp_evec[j][i];
		}
	}

	for (int i = 0; i < n; i++) {

		double norm = 0;

		for (int j = 0; j < n; j++) {

			norm += evec(i, j) * evec(i, j);

		}

		norm = sqrt(norm);

		for (int j = 0; j < n; j++) {
			evec(i, j) = evec(i, j) / norm;
		}
	}

}

void diagonalizeexcmatrix(Matrix<double>& evec, Vektor<double>& eval, Vektor<string>& chlname_vec, Vektor<double>& parameter_vec, Vektor<double> sitevec){

    string str;
    int n = sitevec.size();                         //allocate excitonmatrix that will be diagonalized
    Matrix<double> excitonmatrix(n,n);
    ublas::zero_matrix<double> tmpzeromatrix(n,n);  //to get out chunk
    //cout<<excitonmatrix1<<endl;                   //errorcheking
    excitonmatrix=tmpzeromatrix;


    for(int i=0;i<n;i++){               //set up the diagonal
        excitonmatrix(i,i)=sitevec(i);
    }

    //get the excitonic couplings out from the .dat file
    ifstream exc_coupling ("src/excitonic_coupling.dat");

    double u,v,w;

    if (exc_coupling.is_open()){           //obtain siteenergies from Site_Energies.dat
        while(exc_coupling.good()){           //this is from all other lines
            getline (exc_coupling,str);
            stringstream sstest(str);

            if(str == "NONE"){          //if the user writes MONOMERE into the couplings.dat-file it skips to the Monomere-Option
                break;
            }

            sstest >> u;
            sstest >> v;
            sstest >> w;
            excitonmatrix(u-1,v-1) = w;
            excitonmatrix(v-1,u-1) = w;

        }
    }
    else{   //This part was added so a Monomere can be calculated
            for(int i = 0; i<n; i++){
                for(int j = 0; j<n; j++){
                    if(i!=j){
                        excitonmatrix(i,j) = 0.0;
                    }
                }
            }

    }
//    cout<<"Excitonmatrix Read In: "<<endl;
//    cout<<excitonmatrix<<endl;    //errorchecking

    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Here the Eigenvalueproblem is solved through the NR eigensolver "nr3.h" && eigen_sym
    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //NR MatrixClass .... this is just a temporary object
    MatDoub tmp_excmat(n,n);

    //Copy the excitonmatrix to the NR-Matrixclass

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){
            tmp_excmat[i][j]=excitonmatrix(i,j);
        }
    }

    //NR eigenvalueproblem solver - this creats a class
    Symmeig sol(tmp_excmat);
    MatDoub tmp_evec(n,n);
    VecDoub tmp_eval(n);
	tmp_eval = sol.d;     //obtain solutions from the solverclass
    tmp_evec = sol.z;

    eval.resize(n);         //resize original working objects
    evec.resize(n,n);

    //this is for errorchecking!! the matrix which is diagonalized is just printed out
    /*
    //cout << fixed << setprecision(3);
	cout << "The matrix read in is of size " << n << " x " << n <<" and contains the following elements" << endl << endl;
	// print out the matrix that was read into the program
	for(int i=0;i<n;i++) {
		for(int j=0;j<n;j++) {
			cout << setw(14) << excmat[i][j];
		}
		cout << endl;
	}
	cout << endl;*/


    //Print & Copy all the results
//    cout << endl << "Eigenvalues: " << endl << endl;
    for(int i=0;i<n;i++) {
//		cout << setw(14) << tmp_eval[i];
//		cout << endl;
		eval(i)=tmp_eval[i];
	}

//	cout << endl << "EigenVectors" << endl << endl;
	for(int i=0;i<n;i++) {
//        cout<<"{";
        for(int j=0;j<n;j++){
//            cout << setw(14) << tmp_evec[j][i]<<",";
            evec(i,j)=tmp_evec[j][i];

        }
//        cout<<"}";
//        cout<<endl;
//        cout<<endl;
	}
}

void diagonalize_KineticMatrix(Matrix<double>& A_kinetic, Matrix<double>& evec, Vektor<double>& eval){

    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Here the Eigenvalueproblem is solved through the NR eigensolver "nr3.h" && eigen_unsym
    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //NR MatrixClass .... this is just a temporary object

	int n1 = A_kinetic.size1();
	int n2 = A_kinetic.size2();

	MatDoub tmp_excmat(n1, n2);

    //Copy to the NR-Matrixclass

    for(int i=0; i < n1; i++){
        for(int j=0; j< n2; j++){

        	tmp_excmat[i][j] = A_kinetic(i,j);
        }
    }

    //NR eigenvalueproblem solver - this creates a class
    Unsymmeig sol(tmp_excmat);
    MatDoub tmp_evec(n1, n2);
    VecComplex tmp_eval(n1);
	tmp_eval = sol.wri;     //obtain solutions from the solverclass
    tmp_evec = sol.zz;

    eval.resize(n1);         //resize original working objects
    evec.resize(n1, n2);
    //this is for errorchecking!! the matrix which is diagonalized is just printed out

//    //cout << fixed << setprecision(3);
//	cout << "The matrix read in is of size " << n1 << " x " << n1 <<" and contains the following elements" << endl << endl;
//	// print out the matrix that was read into the program
//	for(int i=0;i<n1;i++) {
//		for(int j=0;j<n1;j++) {
//			cout << setw(14) << tmp_excmat[i][j];
//		}
//		cout << endl;
//	}
//	cout << endl;


    //Print & Copy all the results
	for (int i = 0; i < n1; i++) {
		eval(i) = real(tmp_eval[i]);
	}

	for (int i = 0; i < n1; i++) {
		for (int j = 0; j < n1; j++) {
			evec(i, j) = tmp_evec[j][i];
		}
	}
}

void diagonalize_KineticMatrix_sym(Matrix<double>& A_kinetic, Matrix<double>& evec, Vektor<double>& eval){

    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    //Here the Eigenvalueproblem is solved through the NR eigensolver "nr3.h" && eigen_sym
    //-!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    //NR MatrixClass .... this is just a temporary object
    int n = A_kinetic.size1();
    MatDoub tmp_excmat(n,n);

    //Copy the excitonmatrix to the NR-Matrixclass

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			tmp_excmat[i][j] = A_kinetic(i, j);
		}
	}

    //NR eigenvalueproblem solver - this creats a class
    Symmeig sol(tmp_excmat);
    MatDoub tmp_evec(n,n);
    VecDoub tmp_eval(n);
	tmp_eval = sol.d;     //obtain solutions from the solverclass
    tmp_evec = sol.z;

    eval.resize(n);         //resize original working objects
    evec.resize(n,n);

	//Copy Eigenvalues
	for (int i = 0; i < n; i++) {
		eval(i) = tmp_eval[i];
	}

	//Copy Eigencoefficients
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			evec(i, j) = tmp_evec[j][i];
		}
	}

//	for (int i = 0; i < n; i++) {
//		double norm = 0;
//
//		for (int j = 0; j < n; j++) {
//			norm += evec(i, j) * evec(i, j);
//		}
//
//		norm = sqrt(norm);
//
//		for (int j = 0; j < n; j++) {
//			evec(i, j) = evec(i, j) / norm;
//		}
//	}
}

void get2Exc_Shifts(Matrix<double>& Exc2Shifts, Vektor<double> sitevec){

	int Exc_Mat_Size = sitevec.size();
	Exc2Shifts.resize(Exc_Mat_Size, Exc_Mat_Size);

    for(int i = 0; i < (int)Exc2Shifts.size1(); i++){
    	for(int j = 0; j < (int)Exc2Shifts.size2(); j++){

    		Exc2Shifts(i, j) = 0;
    	}
    }

	ifstream file ("src/Exc2_Shifts.dat");
    string str;

    if (file.is_open()) {
    	while (file.good()) {
    		getline(file, str);

    		if (file) {

    			stringstream Stream(str);

    			int i, j = 0;
    			double Vmn = 0;

    			Stream >> i >> j >> Vmn;

    			Exc2Shifts(i-1, j-1) = Vmn;
    		}
    	}
	}

    for(int i = 0; i < (int)Exc2Shifts.size1(); i++){
    	for(int j = 0; j < (int)Exc2Shifts.size2(); j++){

    		if(i < j){
    			Exc2Shifts(j, i) = Exc2Shifts(i, j);
    		}
    	}
    }
}

#endif // EXCMATRIX_HPP_INCLUDED
