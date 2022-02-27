#ifndef EXCDOMAINS_HPP_INCLUDED
#define EXCDOMAINS_HPP_INCLUDED

#include <vector>
#include <iostream>
#include <iomanip>
#include "Utilities.hpp"
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include "excmatrix.hpp"
#include <boost/assign/std/vector.hpp> // for 'operator+=()'
#include <boost/assert.hpp>;

using namespace boost::assign; // bring 'operator+=()' into scope
using namespace std;

int iselementof(double u, Vektor<double> domainentry){

    int isentry = 0;

    for(int i = 0; i < (int)domainentry.size(); i ++){

            if( u == domainentry(i)){

                isentry = 1;

            }
    }

    return isentry;

}

void deletefromlist(double u, Vektor<double>& domainentry){


    int isentry = 0;
    int position = 0;

    for(int i = 0; i < (int)domainentry.size(); i ++){

            if( u == domainentry(i)){
                isentry = 1;
                position = i;
            }
    }

    if(isentry == 1){
        for(int j = position; j < (int)domainentry.size()-1; j++){

                domainentry(j) = domainentry(j+1);

        }

        domainentry.resize(domainentry.size()-1);
    }

    //else if(isentry == 0){ cout<<"ERROR:void deletefromlist element was not found"<<endl;}

}

void buildexcitondomains_sametype(Vektor<string> chlname_vec){

    //get the excitonic couplings out from the .dat file
    ifstream exc_coupling ("excitonic_coupling_CP29.dat");


    double u,v,w,coupling;
    string str;

    Vektor<double> domCLA(0);
    Vektor<double> domCHL(0);


    for(int j = 0; j < (int)chlname_vec.size(); j++){

        int vecsize;

        if(chlname_vec(j) == "CLA"){

            vecsize = domCLA.size();
            domCLA.resize(vecsize+1);
            domCLA(vecsize) = j+1;

        }

        else if(chlname_vec(j) == "CHL"){

            vecsize = domCHL.size();
            domCHL.resize(vecsize+1);
            domCHL(vecsize) = j+1;

        }

        else{
                cout<<"Error at establishing ExcitonDomains: void buildexcitons()"<<endl;
        }

    }

    cout<<"domCLA: "<<domCLA<<" domCHL: "<<domCHL<<endl;


    Vektor< Vektor<double> > coupling_vec;
    Vektor<double> temp(3);

    if (exc_coupling.is_open()){           //obtain siteenergies from Site_Energies.dat
        while(exc_coupling.good()){           //this is from all other lines
            getline (exc_coupling,str);
            stringstream sstest(str);

            sstest >> u;
            sstest >> v;
            sstest >> coupling;

            temp(0) = u;
            temp(1) = v;
            temp(2) = coupling;


            coupling_vec.resize(coupling_vec.size()+1);
            coupling_vec(coupling_vec.size()-1) = temp;


        }
    }

    //cout<<"Coupling Vector: "<<coupling_vec<<endl;   //errorchecking


    //--------Start with building Exciton domains, seperated into CLA & CHL--------

    Vektor< Vektor<double> > domains;


    //--------Start building CLA Exciton Domains at first Pigment--------

    Vektor<double> CLAavailablepigments;
    CLAavailablepigments.resize(domCLA.size());
    CLAavailablepigments = domCLA;

    cout<<"Available Pigments for CLA Domains: "<<CLAavailablepigments<<endl;

    /*
    deletefromlist(3,availablepigments);    //errorchecking
    cout<<"Available Pigments for Domains: "<<availablepigments<<endl;*/


    while(CLAavailablepigments.size()!=0){

    Vektor<double> domainentry;
    domainentry.resize(domainentry.size()+1);
    domainentry(0) = CLAavailablepigments(0);
    deletefromlist(CLAavailablepigments(0),CLAavailablepigments);

    for(int i = 0; i < (int)domainentry.size(); i++){

        cout<<"i : "<<i<<endl;
        cout<<"Domainentry: "<<domainentry<<endl;

        for(int j = 0; j < (int)coupling_vec.size(); j++){

            if((coupling_vec(j)(0) == domainentry(i))){

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            if(iselementof(u,domCLA)){          //u in excitonic_coupling refers to a CLA

                if(iselementof(v,domCLA)){

                    if(sqrt(coupling*coupling) > 20){       //Coupling > than threshold

                        cout<<u <<" "<<v <<" "<<coupling<<endl;
                        deletefromlist(v,CLAavailablepigments);
                        //cout<<"Available Pigments for Domains: "<<CLAavailablepigments<<endl;

                        if(iselementof(v,domainentry) == 0){
                            domainentry.resize(domainentry.size()+1);
                            domainentry(domainentry.size()-1) = v;
                        }

                    }
                }
            }
            }

            if((coupling_vec(j)(1) == domainentry(i))){

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            if(iselementof(u,domCLA)){          //u in excitonic_coupling refers to a CLA

                if(iselementof(v,domCLA)){

                    if(sqrt(coupling*coupling) > 20){       //Coupling > than threshold

                        cout<<u <<" "<<v <<" "<<coupling<<endl;
                        deletefromlist(u,CLAavailablepigments);
                        //cout<<"Available Pigments for Domains: "<<availablepigments<<endl;

                        if(iselementof(u,domainentry) == 0){
                            domainentry.resize(domainentry.size()+1);
                            domainentry(domainentry.size()-1) = u;
                        }

                    }
                }
            }
            }


        }
    }
    cout<<"Available Pigments for CLA Domains: "<<CLAavailablepigments<<endl;
    cout<<"Domainentry: "<<domainentry<<endl;

    //Sort the Domainentry Vector from Low to High for more beautiful Domain matrices!
    //std::sort is being used

    vector<double> sortvec(domainentry.size());

    for(int z = 0; z < (int)domainentry.size(); z++){
        sortvec[z] = domainentry(z);
        //cout<<"Vectest: "<<sortvec[z]<<endl;
    }

    sort(sortvec.begin(), sortvec.end());

    for(int z = 0; z < (int)domainentry.size(); z++){
        //cout<<"SortedVectest: "<<sortvec[z]<<endl;    //errorchecking
    }

    for(int z = 0; z < (int)domainentry.size(); z++){
        domainentry(z) = sortvec[z];
        //cout<<"Vectest: "<<sortvec[z]<<endl;
    }

    domains.resize(domains.size()+1);
    domains(domains.size()-1) = domainentry;

    cout << domains <<endl;
    }


    //--------Start building CHL Exciton Domains at first Pigment--------

    Vektor<double> CHLavailablepigments;
    CHLavailablepigments.resize(domCHL.size());
    CHLavailablepigments = domCHL;

    cout<<"Available Pigments for CHL Domains: "<<CHLavailablepigments<<endl;

    while(CHLavailablepigments.size()!=0){

    Vektor<double> domainentry;
    domainentry.resize(domainentry.size()+1);
    domainentry(0) = CHLavailablepigments(0);
    deletefromlist(CHLavailablepigments(0),CHLavailablepigments);

    for(int i = 0; i < (int)domainentry.size(); i++){

        cout<<"i : "<<i<<endl;
        cout<<"Domainentry: "<<domainentry<<endl;

        for(int j = 0; j < (int)coupling_vec.size(); j++){

            if((coupling_vec(j)(0) == domainentry(i))){

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            if(iselementof(u,domCHL)){          //u in excitonic_coupling refers to a CLA

                if(iselementof(v,domCHL)){

                    if(sqrt(coupling*coupling) > 20){       //Coupling > than threshold

                        cout<<u <<" "<<v <<" "<<coupling<<endl;
                        deletefromlist(v,CHLavailablepigments);
                        //cout<<"Available Pigments for Domains: "<<CHLavailablepigments<<endl;

                        if(iselementof(v,domainentry) == 0){
                            domainentry.resize(domainentry.size()+1);
                            domainentry(domainentry.size()-1) = v;
                        }

                    }
                }
            }
            }

            if((coupling_vec(j)(1) == domainentry(i))){

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            if(iselementof(u,domCHL)){          //u in excitonic_coupling refers to a CLA

                if(iselementof(v,domCHL)){

                    if(sqrt(coupling*coupling) > 20){       //Coupling > than threshold

                        cout<<u <<" "<<v <<" "<<coupling<<endl;
                        deletefromlist(u,CHLavailablepigments);
                        //cout<<"Available Pigments for Domains: "<<availablepigments<<endl;

                        if(iselementof(u,domainentry) == 0){
                            domainentry.resize(domainentry.size()+1);
                            domainentry(domainentry.size()-1) = u;
                        }

                    }
                }
            }
            }


        }
    }
    cout<<"Available Pigments for CLA Domains: "<<CHLavailablepigments<<endl;
    cout<<"Domainentry: "<<domainentry<<endl;

    //Sort the Domainentry Vector from Low to High for more beautiful Domain matrices!
    //std::sort is being used

    vector<double> sortvec(domainentry.size());
    //Caution: I've used the same sortvector, which was defined above in the CLA part, here for CHL

    for(int z = 0; z < (int)domainentry.size(); z++){
        sortvec[z] = domainentry(z);
        //cout<<"Vectest: "<<vectest[z]<<endl;
    }

    sort(sortvec.begin(), sortvec.end());

    for(int z = 0; z < (int)domainentry.size(); z++){
        //cout<<"SortedVectest: "<<vectest[z]<<endl;    //errorchecking
    }

    for(int z = 0; z < (int)domainentry.size(); z++){
        domainentry(z) = sortvec[z];
        //cout<<"Vectest: "<<sortvec[z]<<endl;
    }

    domains.resize(domains.size()+1);
    domains(domains.size()-1) = domainentry;

    cout << domains <<endl;
    }


    //--------Get the calculated Exciton Domains into matrix form & print them out--------

    cout<<"1st iteration domains obtained: "<<domains<<endl;
    cout<<"There were "<<domains.size()<<" domains buildt"<<endl;

    Vektor< Vektor<double> > finaldomains;

    for(int i = 0; i < (int)domains.size(); i++){

        if(domains(i).size() > 1){

            finaldomains.resize(finaldomains.size()+1);
            finaldomains(finaldomains.size()-1) = domains(i);

        }

    }

    cout<<"Final domains obtained: "<<finaldomains<<endl;
    cout<<"There were "<<finaldomains.size()<<" domains buildt"<<endl;

    Vektor< Matrix<double> > excmatrix_domains;
    excmatrix_domains.resize(finaldomains.size());

    for(int i = 0; i < (int)finaldomains.size(); i++){

        int t = finaldomains(i).size();
        excmatrix_domains(i).resize(t,t);

        for(int j = 0; j < t; j++){

            for(int r = 0; r < t; r++){

                coupling = 0;

                if(j == 0 && r == 0){
                    coupling = 0;
                }

                for(int e = 0; e < (int)coupling_vec.size(); e++){

                    if(coupling_vec(e)(0) == finaldomains(i)(j) && coupling_vec(e)(1) == finaldomains(i)(r)){
                        coupling = coupling_vec(e)(2);
                    }

                    if(coupling_vec(e)(0) == finaldomains(i)(r) && coupling_vec(e)(1) == finaldomains(i)(j)){
                        coupling = coupling_vec(e)(2);
                    }
                }

                excmatrix_domains(i)(j,r) = coupling;
            }
        }
    }

    //Generate File with domains in it for easier overview!

    ofstream domainoutput;
    domainoutput.open ("exciton_domains.out");

    for(int i = 0; i < (int)excmatrix_domains.size(); i++){
        domainoutput<<"Domain "<<i+1<<" generated: "<<endl;

        for(int j = 0; j < (int)excmatrix_domains(i).size1(); j++){
            for(int r = 0; r < (int)excmatrix_domains(i).size2(); r++){

                 domainoutput<<fixed<<setw(15)<<setprecision(5)<<excmatrix_domains(i)(j,r);

            }
            domainoutput<<endl;
        }
    }



    cout<<excmatrix_domains<<endl;

}

void buildexcitondomains(Vektor<string> chlname_vec, Vektor< Matrix<double> >& excmatrix_domains, Vektor< Vektor<double> >& domains, Vektor<double>& parameter_vec){

    //-----Build Exciton-domains-----

    //-----Get the excitonic couplings out from the .dat file-----
    //ifstream exc_coupling ("excitonic_coupling_CP29.dat");
    ifstream exc_coupling ("src/excitonic_coupling.dat");

    double u,v,w,coupling;
    string str;

    Vektor< Vektor<double> > coupling_vec;
    Vektor<double> temp(3);

    if (exc_coupling.is_open()){           //obtain siteenergies from Site_Energies.dat
        while(exc_coupling.good()){           //this is from all other lines
            getline (exc_coupling,str);
            stringstream sstest(str);

            sstest >> u;
            sstest >> v;
            sstest >> coupling;

            temp(0) = u;
            temp(1) = v;
            temp(2) = coupling;


            coupling_vec.resize(coupling_vec.size()+1);
            coupling_vec(coupling_vec.size()-1) = temp;


        }
    }

    //cout<<"Coupling Vector: "<<coupling_vec<<endl;   //errorchecking


    //--------Start with building Exciton domains, no typeseperation is done--------

    //Vektor< Vektor<double> > domains;
    //Vektor< Vektor<double> > domain_siteenergies;

    //--------Start building Exciton Domains at first Pigment--------

    Vektor<double> availablepigments;
    availablepigments.resize(chlname_vec.size());

    for(int i = 0; i < (int)availablepigments.size(); i++){

        availablepigments(i) = i+1;
    }

    //-----Start Iterating-----

    while(availablepigments.size()!=0){                     //until all Pigments from availablepigments are deleted

    Vektor<double> domainentry;
    domainentry.resize(domainentry.size()+1);

    domainentry(0) = availablepigments(0);
    deletefromlist(availablepigments(0), availablepigments);

    for(int i = 0; i < (int)domainentry.size(); i++){            //Start at Pigment 1 (i'th pigment)

//        cout<<"i : "<<i<<endl;
//        cout<<"Domainentry: "<<domainentry<<endl;

        for(int j = 0; j < (int)coupling_vec.size(); j++){       //Go over all other pigments j

            if((coupling_vec(j)(0) == domainentry(i))){     //check coupling vector, built from file

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            //u in excitonic_coupling refers to a CLA

            if(sqrt(coupling*coupling) > parameter_vec(4)){               //Coupling > than threshold

//                cout<<u <<" "<<v <<" "<<coupling<<endl;
                deletefromlist(v,availablepigments);
                        //cout<<"Available Pigments for Domains: "<<CLAavailablepigments<<endl;

                if(iselementof(v,domainentry) == 0){        //Add pigment to domain if it isn't in there already
                    domainentry.resize(domainentry.size()+1);
                    domainentry(domainentry.size()-1) = v;
                }

            }


            }

            if((coupling_vec(j)(1) == domainentry(i))){     //Same loop as above but now couplings from j - i are concerned (because of excitonic_coupling.dat format)

            u = coupling_vec(j)(0);
            v = coupling_vec(j)(1);
            coupling = coupling_vec(j)(2);

            //cout<<"u: "<<u<<" v: "<<v<<" coup: "<<coupling<<endl;

            //u in excitonic_coupling refers to a CLA

            if(sqrt(coupling*coupling) > parameter_vec(4)){       //Coupling > than threshold

//                cout<<u <<" "<<v <<" "<<coupling<<endl;
                deletefromlist(u,availablepigments);
                        //cout<<"Available Pigments for Domains: "<<availablepigments<<endl;

                if(iselementof(u,domainentry) == 0){
                    domainentry.resize(domainentry.size()+1);
                    domainentry(domainentry.size()-1) = u;
                }

            }


            }
        }
    }
//    cout<<"Available Pigments for Domains: "<<availablepigments<<endl;
//    cout<<"Domainentry: "<<domainentry<<endl;

    //Sort the Domainentry Vector from Low to High for more beautiful Domain matrices!
    //std::sort is being used

    vector<double> sortvec(domainentry.size());

    for(int z = 0; z < (int)domainentry.size(); z++){
        sortvec[z] = domainentry(z);
        //cout<<"Vectest: "<<sortvec[z]<<endl;
    }

    sort(sortvec.begin(), sortvec.end());

    for(int z = 0; z < (int)domainentry.size(); z++){
        //cout<<"SortedVectest: "<<sortvec[z]<<endl;    //errorchecking
    }

    for(int z = 0; z < (int)domainentry.size(); z++){
        domainentry(z) = sortvec[z];
        //cout<<"Vectest: "<<sortvec[z]<<endl;
    }

    domains.resize(domains.size()+1);
    domains(domains.size()-1) = domainentry;

    }

    //-----Get the calculated Exciton Domains into matrix form & print them out-----

    cout<<"Domains obtained: "<<domains<<endl;
    cout<<"There were "<<domains.size()<<" domains buildt"<<endl;

    //Vektor< Matrix<double> > excmatrix_domains;
    excmatrix_domains.resize(domains.size());

    for(int i = 0; i < (int)domains.size(); i++){

        int t = domains(i).size();
        excmatrix_domains(i).resize(t,t);

        for(int j = 0; j < t; j++){

            for(int r = 0; r < t; r++){

                coupling = 0;

                if(j == 0 && r == 0){
                    coupling = 0;
                }

                for(int e = 0; e < (int)coupling_vec.size(); e++){

                    if(coupling_vec(e)(0) == domains(i)(j) && coupling_vec(e)(1) == domains(i)(r)){
                        coupling = coupling_vec(e)(2);
                    }

                    if(coupling_vec(e)(0) == domains(i)(r) && coupling_vec(e)(1) == domains(i)(j)){
                        coupling = coupling_vec(e)(2);
                    }
                }

                excmatrix_domains(i)(j,r) = coupling;
            }
        }

    }

//    //-----Correct for the possible P-TrESP Error-----
//    for(int i = 0; i < (int)excmatrix_domains.size(); i++){
//        excmatrix_domains(i) = excmatrix_domains(i) * 1.0;
//    }

    //-----Generate File with domains in it for easier overview!-----

    ofstream domainoutput;
    domainoutput.open ("exciton_domains.out");

	domainoutput << "Domains obtained: " << domains << endl << endl;

	for (int i = 0; i < (int) excmatrix_domains.size(); i++) {
		domainoutput << "Domain " << i + 1 << " generated: " << endl;

		for (int j = 0; j < (int) excmatrix_domains(i).size1(); j++) {
			for (int r = 0; r < (int) excmatrix_domains(i).size2(); r++) {

				domainoutput << fixed << setw(15) << setprecision(5)
						<< excmatrix_domains(i)(j, r);

			}
			domainoutput << endl;
		}
	}

	for (int i = 0; i < (int) excmatrix_domains.size(); i++) {
		cout << "Domain" << i + 1 << " : " << excmatrix_domains(i) << endl;
	}
	domainoutput.close();
}

void build_predefined_excitondomains(Vektor<string> chlname_vec, Vektor< Matrix<double> >& excmatrix_domains, Vektor< Vektor<double> >& domains, Vektor< Vektor<double> >& coupling_vec){

    //-----Build Exciton-domains-----


    //-----Get the excitonic couplings out from the .dat file-----
    //ifstream exc_coupling ("excitonic_coupling_CP29.dat");
    ifstream exc_coupling ("excitonic_coupling.dat");

    double u,v,w,coupling;
    string str;

//    Vektor< Vektor<double> > coupling_vec;
    Vektor<double> temp(3);

    if (exc_coupling.is_open()){           //obtain siteenergies from Site_Energies.dat
        while(exc_coupling.good()){           //this is from all other lines
            getline (exc_coupling,str);
            stringstream sstest(str);

            sstest >> u;
            sstest >> v;
            sstest >> coupling;

            temp(0) = u;
            temp(1) = v;
            temp(2) = coupling;


            coupling_vec.resize(coupling_vec.size()+1);
            coupling_vec(coupling_vec.size()-1) = temp;


        }
    }


//    cout<<"Coupling Vector: "<<coupling_vec<<endl;   //errorchecking
//    exit(1);

//--------Start with building Exciton domains, no typeseperation is done--------

//-----HIER MUSS DIE DOMÄNENGRÖßE GEÄNDERT WERDEN!!-----
//    domains.resize(3);
//    domains.resize(6);
//    vector<int> pre_domain1;
////    pre_domain1 += 1,2,6,7,8,9,10,13;
//
////    pre_domain1 += 1,  4,  16,   19,   22,   25,   28,   37;
//
//    pre_domain1 += 1,  5,   21,   25,   29,   33,   37,   49; //N-Transitions
////    1,  4,  16,   19,   22,   25,   28,   37;
////    1,2,3,  4,5,6,  16,17,18,   19,20,21,   22,23,24,   25,26,27,   28,29,30,   37,38,39;
////    1,2,3,4,   5,6,7,8,   9,10,11,12,     13,14,15,16,    17,18,19,20,    21,22,23,24,    25,26,27,28,    29,30,31,32,    33,34,35,36,    37,38,39,40,    41,42,43,44,     45,46,47,48,   49,50,51,52;
//    vector<int> pre_domain2;
////    pre_domain2 += 3,4,5;
//
////    pre_domain2 += 7,   10,   13
//    pre_domain2 += 9,   13,   17; //N-Transitions
////7,   10,   13;
////    7,8,9,   10,11,12,   13,14,15;
//
//    vector<int> pre_domain3;
////    pre_domain2 += 11,12;
//
////    pre_domain3 += 31,  34;
//    pre_domain3 += 41;   //N-Transitions
////    31,32,33,    34,35,36;
//
//    vector<int> pre_domain4;
//    pre_domain4 += 45;
//
//    vector<int> pre_domain5;
////    pre_domain4 += 2,3, 5,6,    8,9,    11,12,  14,15,  17,18,   20,21,  23,24,   26,27,   29,30,      32,33,      35,36,     38,39,    40,41,42,43,44,45,46,47,48,49,50,51;
//    pre_domain5 += 2,3,4,   6,7,8,   10,11,12,     14,15,16,    18,19,20,    22,23,24,    26,27,28,    30,31,32,    34,35,36,    38,39,40,    42,43,44,     46,47,48,   50,51,52, 53,54,55,56,57,58,59,60,61,62,63,64;
////    pre_domain4 += 2,3, 5,6,    8,9,    11,12,  14,15,  17,18,   20,21,  23,24,   26,27,   29,30,      32,33,      35,36,     38,39,    40,41,42,43,44,45,46,47,48,49,50,51;
//
////    vector<int> pre_domain5;
////    pre_domain5 += 2,3, 5,6,    8,9,    11,12,  14,15,  17,18,   20,21,  23,24,   26,27,   29,30,      32,33,      35,36,     38,39,    40,41,42,43,44,45,46,47,48,49,50,51;
//
////    pre_domain4 += 2,3, 5,6,    8,9,    11,12,  14,15,  17,18,   20,21,  23,24,   26,27,   29,30,      32,33,      35,36,     38,39;
////    pre_domain4 += 40,41,42,43,44,45,46,47,48,49,50,51;
////    pre_domain4 += ;
//
//    vector<int> pre_domain6;
//    pre_domain6 += 1;

    domains.resize(6);
    vector<int> pre_domain1;
    pre_domain1 += 29,2,3,8,9,10,11,12;

    vector<int> pre_domain2;
    pre_domain2 += 4,5,6,7,13,14;

    vector<int> pre_domain3;
    pre_domain3 += 1,16,17,22,23,24,25,26;

    vector<int> pre_domain4;
    pre_domain4 += 18,19,20,21,27,28;;

    vector<int> pre_domain5;
    pre_domain5 += 15,30,31,36,37,38,39,40;

    vector<int> pre_domain6;
    pre_domain6 += 32,33,34,35,41,42;


    domains(0).resize(pre_domain1.size());
    domains(1).resize(pre_domain2.size());
    domains(2).resize(pre_domain3.size());
    domains(3).resize(pre_domain4.size());
    domains(4).resize(pre_domain5.size());
    domains(5).resize(pre_domain6.size());

    for(int i = 0; i < (int)pre_domain1.size(); i++){
            //cout<< "values : "<< pre_domain1[i] <<endl;
            domains(0)(i) = pre_domain1[i];
    }
    for(int i = 0; i < (int)pre_domain2.size(); i++){
//            cout<< "values : "<< pre_domain2[i] <<endl;
            domains(1)(i) = pre_domain2[i];
    }
    for(int i = 0; i < (int)pre_domain3.size(); i++){
//            cout<< "values : "<< pre_domain3[i] <<endl;
            domains(2)(i) = pre_domain3[i];
    }

    for(int i = 0; i < (int)pre_domain4.size(); i++){
//            cout<< "values : "<< pre_domain4[i] <<endl;
            domains(3)(i) = pre_domain4[i];
    }

    for(int i = 0; i < (int)pre_domain5.size(); i++){
//            cout<< "values : "<< pre_domain4[i] <<endl;
            domains(4)(i) = pre_domain5[i];
    }

    for(int i = 0; i < (int)pre_domain6.size(); i++){
//            cout<< "values : "<< pre_domain4[i] <<endl;
            domains(5)(i) = pre_domain6[i];
    }


    excmatrix_domains.resize(domains.size());

    for(int i = 0; i < (int)domains.size(); i++){

        int t = domains(i).size();
        excmatrix_domains(i).resize(t,t);

        for(int j = 0; j < t; j++){

            for(int r = 0; r < t; r++){

                coupling = 0;

                if(j == 0 && r == 0){
                    coupling = 0;
                }

                for(int e = 0; e < (int)coupling_vec.size(); e++){

                    if(coupling_vec(e)(0) == domains(i)(j) && coupling_vec(e)(1) == domains(i)(r)){
                        coupling = coupling_vec(e)(2);
                    }

                    if(coupling_vec(e)(0) == domains(i)(r) && coupling_vec(e)(1) == domains(i)(j)){
                        coupling = coupling_vec(e)(2);
                    }
                }

                excmatrix_domains(i)(j,r) = coupling;
            }
        }

    }



    //Correct for the possible P-TrESP Error
    for(int i = 0; i < (int)excmatrix_domains.size(); i++){

        excmatrix_domains(i) = excmatrix_domains(i) * 1.0;
    }


    //Generate File with domains in it for easier overview!

    ofstream domainoutput;
    domainoutput.open ("exciton_domains.out");

    domainoutput<<"Domains obtained: "<<domains<<endl<<endl;

    for(int i = 0; i < (int)excmatrix_domains.size(); i++){
        domainoutput<<"Domain "<<i+1<<" generated: "<<endl;

        for(int j = 0; j < (int)excmatrix_domains(i).size1(); j++){
            for(int r = 0; r < (int)excmatrix_domains(i).size2(); r++){

                 domainoutput<<fixed<<setw(15)<<setprecision(5)<<excmatrix_domains(i)(j,r);

            }
            domainoutput<<endl;
        }
    }

    for(int i = 0; i < (int)excmatrix_domains.size(); i++){
	cout<<"Domain"<<i+1<<" : "<<excmatrix_domains(i)<<endl;
    }
    domainoutput.close();
}

void addsiteenergies(Vektor< Matrix<double> >& excmatrix_domains, Vektor< Vektor<double> >& domains, Vektor<double>& sitevec){

    //cout<<"Domains obtained: "<<domains<<endl;
    //cout<<"There were "<<domains.size()<<" domains buildt"<<endl;

    for(int i = 0; i < (int)domains.size(); i++){

        for(int j = 0; j < (int)domains(i).size(); j++){

            excmatrix_domains(i)(j,j) = sitevec(domains(i)(j)-1);
        }
    }
    //cout<<"Excmatrix after adding Siteenergies: "<<excmatrix_domains<<endl;       //errorchecking

}

Vektor< Matrix<double> > addsiteenergies_parallel(Vektor< Matrix<double> > excmatrix_domains, Vektor< Vektor<double> > domains, Vektor<double> sitevec){

    //cout<<"Domains obtained: "<<domains<<endl;
    //cout<<"There were "<<domains.size()<<" domains buildt"<<endl;

    Vektor< Matrix<double> > temp_excmatrix_domains;
    temp_excmatrix_domains = excmatrix_domains;

    for(int i = 0; i < (int)domains.size(); i++){

        for(int j = 0; j < (int)domains(i).size(); j++){

            temp_excmatrix_domains(i)(j,j) = sitevec(domains(i)(j)-1);
        }
    }

    return temp_excmatrix_domains;
    //cout<<"Excmatrix after adding Siteenergies: "<<excmatrix_domains<<endl;       //errorchecking

}

double getcouplings(Vektor< Vektor<double> >& coupling_vec, int u, int v){

    for(unsigned long i = 0; i < coupling_vec.size(); i++){
        if((coupling_vec(i)(0) == u) && (coupling_vec(i)(1) == v)){
            return coupling_vec(i)(2);
        }
    }

    for(unsigned long i = 0; i < coupling_vec.size(); i++){
        if((coupling_vec(i)(0) == v) && (coupling_vec(i)(1) == u)){
            return coupling_vec(i)(2);
        }
    }

    cout << "NO SUITING COUPLING FOUND IN GETCOUPLINGS FUNCTION!" <<endl;
    return -9999999999999999;
}

void stoerungstheorie(Matrix<double>& evecs, Vektor<double>& evals, Vektor<double> evalshigher, Matrix<double> evecshigher, Matrix< double > coupling_mat, Vektor< Vektor<double> > domains, int domainnumber, int domainstoerung){

    int evalssize = evals.size();
    int evecssize = evecs.size1();
    int evecssize2 = evecs.size2();
    int evecshighersize = evecshigher.size1();


    //*****Evals Erweiterung*****

    evals.resize((int)evals.size() + (int)evalshigher.size());

    for(int i = 0; i < (int)evalshigher.size(); i++){

        evals(i + evalssize) = evalshigher(i);
    }

    //*****Die Excitonenzustände M für Qy*****

    evecs.resize(evecs.size1(), evecs.size2() + (int)evalshigher.size());

    for(int M = 0; M < (int)evalssize; M++){

        for(int n = 0; n < (int)evalshigher.size(); n++){


            double V = 0.0;

            for(int N = 0; N < (int)evalshigher.size(); N++){

                double cn = 0.0;

                for(int i = 0; i < (int)evecs.size1(); i++){               //counts c(M)

                    for(int j = 0; j < (int)evecshigher.size1(); j++){     //counts c(N)

                        cn += evecs(M,i) * evecshigher(N,j) * coupling_mat(domains(domainnumber)(i),domains(domainstoerung)(j));

//                        cout << i <<"  "<< j << "  " << evecs(M,i) << "  " << evecshigher(N,j) << "  " << getcouplings(coupling_vec,domains(domainnumber)(i),domains(domainstoerung)(j)) <<endl;
                    }
                }

                cn = cn / (evals(M) - evalshigher(N)) * evecshigher(N,n);
//                cout << "n = " << n << "  cn = " << cn << " N = " << N <<"  M = " << M <<endl;

                V += cn;
            }

            evecs(M,evalssize + n) = V;
        }
    }

//    cout << "Evecs in function Störungstheorie : " <<endl;
//    cout << evecs <<endl;


    //*****Die Excitonenzustände dM für höhere Zustände*****

    evecs.resize(evecs.size1() + evalshigher.size(), evecs.size2());


    for(int N = 0; N < (int)evalshigher.size(); N++){

        for(int m = 0; m < (int)evalssize; m++){


            double V = 0.0;

            for(int M = 0; M < (int)evalssize; M++){

                double cn = 0.0;

                for(int i = 0; i < (int)evecshigher.size1(); i++){               //counts c(M)

                    for(int j = 0; j < (int)evecssize; j++){     //counts c(N)

                       cn += evecshigher(N,i) * evecs(M,j) * coupling_mat(domains(domainnumber)(j),domains(domainstoerung)(i));

                    }
                }

                cn = cn / (evalshigher(N) - evals(M)) * evecs(M,m);
//                cout << "n = " << n << "  cn = " << cn << " N = " << N <<"  M = " << M <<endl;

                V += cn;
            }

            evecs(evalssize + N,m) = V;
        }
    }

//    cout <<" Evecs : " << endl;
//    cout << evecs <<endl;

    //*****Die Excitonenzustände dM für höhere Zustände auffüllen mit ihren Eigencoeffizienten*****

    for(int M = 0; M < evecshighersize; M++){
        for(int N = 0; N < evecshighersize; N++){

            evecs(evecssize + M, evecssize + N) = evecshigher(M,N);
        }
    }

    for(int i = 0; i < (int)evecs.size1(); i++){
        double norm = 0;

        for(int j = 0; j < (int)evecs.size1(); j++){

            norm += evecs(i,j) * evecs(i,j);

        }

        norm = sqrt(norm);

        for(int j = 0; j < (int)evecs.size1(); j++){
            evecs(i,j) = evecs(i,j) / norm;
        }
    }

//    cout <<" Evecs Final : " << endl;
//    cout << evecs <<endl;
//
//    cout <<" Evals Final : " << endl;
//    cout << evals <<endl;


//    cout<<"Domains: "<<endl;
//    cout<<"Getting Couplings"<<endl;
//
//    cout << domains(domainnumber)(0) << "    " << domains(domainstoerung)(2)<<endl;
//    cout << getcouplings(coupling_vec, domains(domainnumber)(0),domains(domainstoerung)(2))<<endl;
////    cout << getcouplings(coupling_vec, domains(domainnumber)(0),domains(domainstoerung)(0))<<endl;
//    cout<<domains(0)(1)<<endl;

}

void stoerungstheorie_singlestate(Matrix<double>& evecs, Vektor<double>& evals, Vektor<double> evalshigher, Vektor< Vektor<double> > coupling_vec, Vektor< Vektor<double> > domains, int domainnumber, int domainstoerung){
	//*****Störungstheorie für einzelne Pigmente*****

    int evalssize = evals.size();
    int evecssize = evecs.size1();
//    int evecssize2 = evecs.size2();

//    cout << evecssize << " size 2: " << evecssize2 <<endl;

    evals.resize((int)evals.size() + (int)evalshigher.size());

    for(int i = 0; i < (int)evalshigher.size(); i++){

        evals(i + evalssize) = evalshigher(i);
    }

    cout << "Evals new: "<<endl<<endl;
    cout << evals <<endl;

    //*****Die Excitonenzustände M für Qy*****

    evecs.resize(evecs.size1(), evecs.size2() + evalshigher.size());

    for(int M = 0; M < evecssize; M++){

        for(int j = 0; j < (int)evalshigher.size(); j++){

            double Vab = 0;

            for(int i = 0; i < evecssize; i++){

                    Vab += evecs(M,i) * getcouplings(coupling_vec,domains(domainnumber)(i),domains(domainstoerung)(j));
            }

            Vab = Vab / (evals(M) - evalshigher(j));

            evecs(M,evecssize + j) = Vab;
        }
    }

    //*****Die Excitonenzustände dM für höhere Zustände*****

    evecs.resize(evecs.size1() + evalshigher.size(), evecs.size2());

    for(int j = 0; j < (int)evalshigher.size(); j++){

        for(int M = 0; M < evecssize; M++){

                double Vab = 0;

                for(int i = 0; i < evecssize; i++){

                    Vab += evecs(M,i) * getcouplings(coupling_vec,domains(domainnumber)(i),domains(domainstoerung)(j));
                }

                Vab = Vab / (evalshigher(j) - evals(M));
//        cout << Vab <<endl;

                evecs(evecssize + j,M) = Vab;
        }
    }

    //*****Die erweiterung von evecs auf 0.0 setzen*****

    for(int i = 0; i < (int)evalshigher.size();i++){
        for(int j = 0; j < (int)evalshigher.size();j++){

            evecs(evecssize + i, evecssize + j) = 0.0;
        }
    }

    //    cout <<" Evecs : " << endl;
    //    cout << evecs <<endl;

    //*****Die Excitonenzustände dM für höhere Zustände auffüllen mit ihren Eigencoeffizienten*****

    for(int j = 0; j < (int)evalshigher.size(); j++){

        evecs(evecssize + j, evecssize + j) = 1.0;
    }


    for(int i = 0; i < (int)evecs.size1(); i++){
        double norm = 0;

        for(int j = 0; j < (int)evecs.size1(); j++){

            norm += evecs(i,j) * evecs(i,j);

        }

        norm = sqrt(norm);

        for(int j = 0; j < (int)evecs.size1(); j++){
            evecs(i,j) = evecs(i,j) / norm;
        }
    }


//    cout <<" Evecs Final : " << endl;
//    cout << evecs <<endl;
//
//    cout <<" Evals Final : " << endl;
//    cout << evals <<endl;


//    cout<<"Domains: "<<endl;
//    cout<<"Getting Couplings"<<endl;
//
//    cout << domains(domainnumber)(0) << "    " << domains(domainstoerung)(2)<<endl;
//    cout << getcouplings(coupling_vec, domains(domainnumber)(0),domains(domainstoerung)(2))<<endl;
////    cout << getcouplings(coupling_vec, domains(domainnumber)(0),domains(domainstoerung)(0))<<endl;
//    cout<<domains(0)(1)<<endl;

}

void ExcitonTransitionDipoles(Vektor< Vektor<double> >& trans_dipole_moment, Vektor< Vektor<double> >& trans_dipole_vector, Vektor<double> evalshigher, Matrix<double> evecshigher, Vektor< Vektor< Vektor<double> > > dipolemat){

    //ABSORPTION-PART
//    Vektor< Vektor<double> > trans_dipole_moment(evalshigher.size());

    trans_dipole_moment.resize(evalshigher.size());
    trans_dipole_vector.resize(evalshigher.size());

    ublas::zero_vector<double> zero(3);
    for(int i=0;i<(int)trans_dipole_moment.size();i++){

        trans_dipole_moment(i) = zero;    //zerovectors so += can be used!
        trans_dipole_vector(i) = zero;
    }

    for(int i = 0; i<(int)trans_dipole_moment.size(); i++){              //check if error occurs!!
        for(int j = 0; j<(int)trans_dipole_moment.size(); j++){

            trans_dipole_moment(i) += dipolemat(0)(j) * evecshigher(i,j);
            trans_dipole_vector(i) += dipolemat(1)(j) * evecshigher(i,j);       //careful! dipolemat(0) because it contains the ei_dipole
        }                                                                       //there are only 2 entries until now!!!
    }

}

void MakeCouplingMatrix(Vektor< Vektor<double> > coupling_vec, Matrix<double>& coupling_mat, unsigned long matrixsize){

    coupling_mat.resize(matrixsize, matrixsize);

    for(unsigned long i = 0; i < coupling_mat.size1(); i++){
        for(unsigned long j = 0; j < coupling_mat.size2(); j++){
            coupling_mat(i,j) = 0.0;
        }
    }

    for(unsigned long i = 0; i < coupling_vec.size(); i++){

        coupling_mat(coupling_vec(i)(0), coupling_vec(i)(1)) = coupling_vec(i)(2);
    }

    for(unsigned long i = 0; i < coupling_mat.size1(); i++){
        for(unsigned long j = 0; j < coupling_mat.size2(); j++){
            coupling_mat(j,i) = coupling_mat(i,j);
        }
    }

//    ofstream test;
//    test.open("makecouplingmat.txt");
//    test << coupling_mat;

}

void Build_2Exiton_Matrix(Matrix<double>& Exciton_Matrix, Matrix<double>&  Exciton2_Matrix){

	int Exc_Mat_Size = Exciton_Matrix.size1();
	int Exc2_Mat_Size = Exc_Mat_Size * (Exc_Mat_Size - 1) / 2;

	Exciton2_Matrix.resize(Exc2_Mat_Size, Exc2_Mat_Size);

	for(int s1 = 0; s1 < (int)Exciton2_Matrix.size1(); s1++){
		for(int s2 = 0; s2 < (int)Exciton2_Matrix.size2(); s2++){
			Exciton2_Matrix(s1, s2) = 0;
		}
	}

	//-----2-Exciton Matrix Diagonale-----
	int k = 0;

	for (int i = 0; i < Exc_Mat_Size - 1; i++) {
		for (int j = i + 1; j < Exc_Mat_Size; j++) {

			Exciton2_Matrix(k, k) = Exciton_Matrix(i,i) + Exciton_Matrix(j,j);
			k++;
//			cout << "i: " << i << "  j: " << j << endl;
		}
	}

//	Exc_Mat_Size = 3;

//	cout << "Building 2-Exciton Off-Diagonal" << endl;

	int u = 0;
	int v = 0;


	for (int i = 0; i < Exc_Mat_Size - 1; i++) {
		for (int j = i + 1; j < Exc_Mat_Size; j++){

			for(int m = 0; m < Exc_Mat_Size - 1; m++){
				for(int n = m + 1; n < Exc_Mat_Size; n++){


					if(m == i && n != j){

						Exciton2_Matrix(u,v) = Exciton_Matrix(n,j);
					}
					else if(m == j && n != i){

						Exciton2_Matrix(u,v) = Exciton_Matrix(n,i);
					}
					else if(n == i && m != j){

						Exciton2_Matrix(u,v) = Exciton_Matrix(m,j);
					}
					else if(n == j && m != i){

						Exciton2_Matrix(u,v) = Exciton_Matrix(m,i);
					}

					u++;

//					cout << "All: i j m n: " << i+1 << " " << j+1 << " " << m+1 << " "  << n+1 << endl;
//
//					if(m == i && n != j){cout << " j:" << j+1 << " n: " << n+1 << endl;}
//
//					if(m == j && n != i){cout << " i:" << i+1 << " n: " << n+1 << endl;}
//
//					if(n == i && m != j){cout << " j:" << j+1 << " m: " << m+1 << endl;}
//
//					if(n == j && m != i){cout << " i:" << i+1 << " m: " << m+1 << endl;}

				}
			}
			u = 0;
			v++;
		}
	}
}

void Build_2Exiton_Matrix_DEBUG(Matrix<double>& Exciton_Matrix, Matrix<double>&  Exciton2_Matrix, Matrix<double>& Exc2Shifts){

	int Exc_Mat_Size = Exciton_Matrix.size1();
	int Exc2_Mat_Size = Exc_Mat_Size * (Exc_Mat_Size - 1) / 2;

	Exciton2_Matrix.resize(Exc2_Mat_Size, Exc2_Mat_Size);

	for(int s1 = 0; s1 < (int)Exciton2_Matrix.size1(); s1++){
		for(int s2 = 0; s2 < (int)Exciton2_Matrix.size2(); s2++){
			Exciton2_Matrix(s1, s2) = 0;
		}
	}

	//-----2-Exciton Matrix Diagonale-----
	int k = 0;

	for (int i = 0; i < Exc_Mat_Size - 1; i++) {
		for (int j = i + 1; j < Exc_Mat_Size; j++) {

			Exciton2_Matrix(k, k) = Exciton_Matrix(i,i) + Exciton_Matrix(j,j) + Exc2Shifts(i, j);
			k++;
		}
	}

}

void Build_2Exiton_Matrix_Shifts(Matrix<double>& Exciton_Matrix, Matrix<double>&  Exciton2_Matrix, Matrix<double>& Exc2Shifts){

	int Exc_Mat_Size = Exciton_Matrix.size1();
	int Exc2_Mat_Size = Exc_Mat_Size * (Exc_Mat_Size - 1) / 2;

	Exciton2_Matrix.resize(Exc2_Mat_Size, Exc2_Mat_Size);

	for(int s1 = 0; s1 < (int)Exciton2_Matrix.size1(); s1++){
		for(int s2 = 0; s2 < (int)Exciton2_Matrix.size2(); s2++){
			Exciton2_Matrix(s1, s2) = 0;
		}
	}

	//-----2-Exciton Matrix Diagonale-----
	int k = 0;

	for (int i = 0; i < Exc_Mat_Size - 1; i++) {
		for (int j = i + 1; j < Exc_Mat_Size; j++) {

			Exciton2_Matrix(k, k) = Exciton_Matrix(i,i) + Exciton_Matrix(j,j) + Exc2Shifts(i, j);
			k++;
		}
	}

	int u = 0;
	int v = 0;


	for (int i = 0; i < Exc_Mat_Size - 1; i++) {
		for (int j = i + 1; j < Exc_Mat_Size; j++){

			for(int m = 0; m < Exc_Mat_Size - 1; m++){
				for(int n = m + 1; n < Exc_Mat_Size; n++){


					if(m == i && n != j){

						Exciton2_Matrix(u,v) = Exciton_Matrix(n,j);
					}
					else if(m == j && n != i){

						Exciton2_Matrix(u,v) = Exciton_Matrix(n,i);
					}
					else if(n == i && m != j){

						Exciton2_Matrix(u,v) = Exciton_Matrix(m,j);
					}
					else if(n == j && m != i){

						Exciton2_Matrix(u,v) = Exciton_Matrix(m,i);
					}

					u++;
				}
			}
			u = 0;
			v++;
		}
	}
}

void Build_2Exc_Evec_Map(Matrix<double> evecs, Matrix<double>& Exc2_Mat){

	int Evecs_Size = evecs.size1();
//	int Evecs_Size = 5;
	int Exc2_Mat_Size = Evecs_Size * (Evecs_Size - 1) / 2;

	Exc2_Mat.resize(Evecs_Size, Evecs_Size);

	for(int i = 0; i < (int)Exc2_Mat.size1(); i++){
		for(int j = 0; j < (int)Exc2_Mat.size2(); j++){
			Exc2_Mat(i,j) = 0.0;
		}
	}

	int k = 0;
//	cout << "Mapping Matrix from 2-Exciton State |mn> to 1D: {l1,l2,l3,...}" << endl;

	for (int i = 0; i < Evecs_Size - 1; i++) {
		for (int j = i + 1; j < Evecs_Size; j++) {

			Exc2_Mat(i,j) = k;
			k++;
//			cout << "i: " << i << "  j: " << j << endl;
		}
	}

	for (int i = 0; i < Evecs_Size - 1; i++) {
		for (int j = i + 1; j < Evecs_Size; j++) {

			Exc2_Mat(j, i) = Exc2_Mat(i, j);
		}
	}

//	cout << "Exc2_Mat: " << Exc2_Mat << endl;

}

double Foerster_Couplings(int M_D1, int N_D2, Matrix<double>& evecs_D1, Matrix<double>& evecs_D2, Vektor< Vektor<double> >& coupling_vec, Vektor< Vektor<double> >& domains, int dom_a, int dom_b){

	double sum = 0;


	for (int m = 0; m < (int) domains(dom_a).size(); m++) {
		for (int n = 0; n < (int) domains(dom_b).size(); n++) {

			sum += evecs_D1(M_D1, m) * evecs_D2(N_D2, n) * getcouplings(coupling_vec, domains(dom_a)(m), domains(dom_b)(n));

//			cout << "Foerster Coupling: " << m << "," << n << "  " << evecs_D1(M_D1, m) << "  " << evecs_D2(N_D2, n) << "  " << getcouplings(coupling_vec, domains(dom_a)(m), domains(dom_b)(n)) << "  " << sum << endl;
//			cout << "Domains: " << domains(dom_a)(m) << "  " << domains(dom_b)(n) << endl;
		}
	}

	return sum;

}

void Get_Coupling_Vec(Vektor< Vektor<double> >& coupling_vec){

	//-----Get the full excitonic couplings out from the .dat file-----
    ifstream exc_coupling ("src/excitonic_coupling.dat");

    double u,v,w,coupling;
    string str;

    Vektor<double> temp(3);

    if (exc_coupling.is_open()){
        while(exc_coupling.good()){
            getline (exc_coupling,str);
				if (exc_coupling){
				stringstream sstest(str);

				sstest >> u;
				sstest >> v;
				sstest >> coupling;

				temp(0) = u;
				temp(1) = v;
				temp(2) = coupling;


				coupling_vec.resize(coupling_vec.size()+1);
				coupling_vec(coupling_vec.size()-1) = temp;
            }

        }
    }


//    cout<<"Coupling Vector: "<<coupling_vec<<endl;   //errorchecking
}

void Build_Huang_Domains(Vektor<double>& Huang_Fac, Vektor< Vektor<double> >& domains, Vektor< Vektor<double> >& Huang_Domains){

	//-----Build a Matrix similar to domains but with Huang-Factors of the pigments instead of just the pigment number (as in domains)-----

	Huang_Domains = domains;

	for(int i = 0; i < (int)domains.size(); i++){
		for(int j = 0; j < (int)domains(i).size();j++){

			int Pig_Numb = domains(i)(j);
			Huang_Domains(i)(j) = Huang_Fac(Pig_Numb - 1);
//			Huang_Domains(i)(j) = 1.0;
		}
	}

}

#endif // EXCDOMAINS_HPP_INCLUDED
