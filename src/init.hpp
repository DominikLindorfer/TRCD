#ifndef INIT_HPP_INCLUDED
#define INIT_HPP_INCLUDED

#include <bits/stdc++.h>
#include "Utilities.hpp"
#include "ChlClass.hpp"

using namespace std;

void createNBNDfile(char *pdbname){

    //Open .pdb file and obtain all NB-ND coordinates
    ifstream pdbfile (pdbname);

    //allign a NB-ND file with all the corresponding coordinates
    ofstream nbndfile;
    nbndfile.open("nb_nd.dat",ios::out);



    //create loopstring & strings for stringsearch
    string str;
    string str2 ("NB  C");
    string str3 ("ND  C");
    string temp;

    if (pdbfile.is_open()){
        while (pdbfile.good()){
            getline (pdbfile,str);  //gets every line until EOF

        if (str.find(str2) != string::npos){
            //cout<<str<<endl;

            unsigned pos = str.find("1.00");    //maybe needs to be changed for other files
            unsigned pos1 = str.find("NB");

            temp=str.substr (pos1,pos-pos1);    //create substring
            //cout<<temp<<endl;
            nbndfile<<temp<<endl;
            }

        if (str.find(str3) != string::npos){
            //cout<<str<<endl;

            unsigned pos = str.find("1.00");     //maybe needs to be changed for other files
            unsigned pos1 = str.find("ND");

            temp=str.substr (pos1,pos-pos1);    //create substring
            //cout<<temp<<endl;
            nbndfile<<temp<<endl;
            }
        }
        pdbfile.close();
    }

    else cout << "Unable to open .pdb file";
}

void init(Vektor<double>& debye_vec,Vektor<double>& parameter_vec, Vektor<double>& Huang_Fac, Vektor<double>& Inh_Broads){

    //Init.dat will be analyzed
    ifstream init ("src/init.dat");
    string str;
    string pdbfilename;
    string tmp;
    double temp;

    //First Line of init.dat is pdbfile name!
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word
            sstest >> pdbfilename;
        }
    }

    //Second Line checks if NB-ND file needs to be set up
    if (init.is_open()){
        if (init.good()){
            getline (init,str);

            if (str.compare(0,1,"1") == 0){        //check if filecreation is enabled (1) or not (0)
                cout<<"NB-ND File will be created!"<<endl;
                cout<<".pdb filename is: "<<pdbfilename<<endl;

                char buffer[200];
                size_t length = pdbfilename.copy(buffer,str.length(),0);
                buffer[length]='\0';
                createNBNDfile(buffer);
            }
        }
    }

    //Third Line checks Dipolestrength of CHL & CLA     //this might need some changes to fit the dipolematrix
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word

            int k=0;                    //get every
            while(sstest >> temp){

                debye_vec.resize(debye_vec.size()+1);
                debye_vec(k)=temp;
                k++;
            }
            //cout<<debye_vec<<endl;    //errorchecking
        }
    }
    //Fourth Line gets set Temperature into the Parameter-Vector, parameter_vec
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word
            parameter_vec.resize(parameter_vec.size()+1);
            sstest >> parameter_vec(0);
            //cout<<debye_vec<<endl;    //errorchecking
        }
    }
    //Fifth Line gets the Vibration-Correlation Radius Rc into the Parameter-Vector, parameter_vec
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word
            parameter_vec.resize(parameter_vec.size()+1);
            sstest >> parameter_vec(1);
            //cout<<parameter_vec<<endl;    //errorchecking
        }
    }
    //Sixth Line get the FWHM for the static disorder
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word
            parameter_vec.resize(parameter_vec.size()+1);
            sstest >> parameter_vec(2);
            //cout<<parameter_vec<<endl;    //errorchecking
        }
    }

    //Seventh Line gets the number of static disorder loops
    if (init.is_open()){
        if (init.good()){
            getline (init,str);
            stringstream sstest(str);   //convert string to stream to get out only first word
            parameter_vec.resize(parameter_vec.size()+1);
            sstest >> parameter_vec(3);
            //cout<<parameter_vec<<endl;    //errorchecking
        }
    }

    //Eighth Line gets the Exciton Domain Radius
    if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(4);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

    //Ninth Line gets the Pump Pulse Width
    if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(5);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

    //Tenth Line gets the Pump Pulse Wavelength
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(6);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

	//Eleventh Line gets the Delay-Time between Pump and Probe Pulse
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(7);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

	//Twelfth Line gets the Huang-Rhys Factors
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);

			stringstream sstest(str);   //convert string to stream to get out only first word

			int Diff_Huang = 0;
			sstest >> Diff_Huang;

			for(int i = 0; i < Diff_Huang; i++){

				int Pig_Numb = 0;
				int old_size = Huang_Fac.size();
				double Huang_Value = 0;

				sstest >> Pig_Numb;
				sstest >> Huang_Value;

				Huang_Fac.resize(old_size + Pig_Numb);

				for(int j = old_size; j < (int)Huang_Fac.size(); j++){
					Huang_Fac(j) = Huang_Value;
				}

			}
		}
	}

	//Thirteenth Line gets the Inhomogeneous Broadenings
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);

			stringstream sstest(str);

			Vektor<double> Inh_Broads_Diff;

			while (sstest.good()){
				double value;
				sstest >> value;
				push_back(Inh_Broads_Diff, value);
			}
			//-----Delete the last entry which should be 0 -----
			Inh_Broads_Diff.resize(Inh_Broads_Diff.size() - 1);

			//-----Set the Inh_Broads Vektor-----
			Inh_Broads.resize(Huang_Fac.size());
			for(int i = 0; i < (int)Inh_Broads.size(); i++){
				Inh_Broads(i) = parameter_vec(2);
			}
			for(int i = 0; i < (int)Inh_Broads_Diff.size(); i += 2){
				Inh_Broads( Inh_Broads_Diff(i) - 1 ) = Inh_Broads_Diff(i + 1);
			}
		}
	}


    //Fourteenth line gets the Cutoff for Kinetic Matrix
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(8);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

    //Fifteenth line gets the CLA Fluroescence Time
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(9);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

    //Sixteenth line gets the RC e- Transferrate / Time
	if (init.is_open()) {
		if (init.good()) {
			getline(init, str);
			stringstream sstest(str);   //convert string to stream to get out only first word
			parameter_vec.resize(parameter_vec.size() + 1);
			sstest >> parameter_vec(10);
			//cout<<parameter_vec<<endl;    //errorchecking
		}
	}

}

void getdipoles(Vektor< Vektor< Vektor<double> > >& dipolemat, Vektor<double> debye_vec, Vektor<string>& chlname_vec){

    //At first let us obtain the unitvectors of the dipole transitionmoments from the NB-ND file
    //they will be stored as vectors within an Vektor

    Vektor< Vektor<double> > ei_dipole(0);    		//Dipolevector of different ei-unitvectors
    Vektor< Vektor<double> > temp_ei_dipole;
    Vektor< Vektor<double> > RCenter_dipole(0);  	//Center of Dipole

    ifstream nbndfile ("src/nb_nd.dat");
    string str; //dummystring
    string tmpstr;


    int k=0;
    Vektor<double> ri_dipole(3);
    Vektor<double> rj_dipole(3);

	if (nbndfile.is_open()) {
		while (nbndfile.good()) {
			getline(nbndfile, str);
			if (nbndfile) {
				stringstream sstest(str);

				sstest >> tmpstr;

				chlname_vec.resize(chlname_vec.size() + 1);
				sstest >> chlname_vec(k);                   //get CHL || CLA || w/e into the namevector

				for (int i = 2; i < 4; i++) {
					sstest >> tmpstr;
				}
				for (int i = 0; i < 3; i++) {
					sstest >> ri_dipole(i);
				}

				getline(nbndfile, str);
				stringstream sstest1(str);

				for (int i = 0; i < 4; i++) {
					sstest1 >> tmpstr;
				}
				for (int i = 0; i < 3; i++) {
					sstest1 >> rj_dipole(i);
				}

				ei_dipole.resize(ei_dipole.size() + 1);
				ei_dipole(k) = (rj_dipole - ri_dipole) / norm_2(rj_dipole - ri_dipole);

				RCenter_dipole.resize(RCenter_dipole.size() + 1);
				RCenter_dipole(k) = 1.0 / 2.0 * (ri_dipole + rj_dipole);

				k++;
			}
			else
				break;
		}
	}
    //cout<<ei_dipole<<endl;          //errorchecking
    //cout<<RCenter_dipole<<endl;     //errorchecking

	//Dipolematrix has UnitTransition Dipolevectors in & The Centers of the CHL for RM
	dipolemat.resize(dipolemat.size() + 2);

	dipolemat(0) = ei_dipole;
	dipolemat(1) = RCenter_dipole;

    //sets the dipolemoments according to init.dat debye_a, debye_b

    //cout<<"Dipolemat size: "<<dipolemat(0).size()<<endl;      //errorchecking
    //cout<<"Debye_Vec size: "<<debye_vec.size()<<endl;
    //cout<<"Chl_Namevec size: "<<chlname_vec.size()<<endl;
    //cout<<chlname_vec<<endl;
    //cout<<dipolemat(0)<<endl;

    for(int i = 0; i < (int)dipolemat(0).size(); i++){
        if(chlname_vec(i) == "CLA"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(0);
        }
        if(chlname_vec(i) == "CHL"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(1);
        }
        if(chlname_vec(i) == "CLABy"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(2);
        }
        if(chlname_vec(i) == "CHLBy"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(3);
        }
        if(chlname_vec(i) == "CLAQx"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(4);
        }
        if(chlname_vec(i) == "CHLQx"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(5);
        }
        if(chlname_vec(i) == "CLABx"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(6);
        }
        if(chlname_vec(i) == "CHLBx"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(7);
        }
        if(chlname_vec(i) == "LUT"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(8);
        }
        if(chlname_vec(i) == "LUT01"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(9);
        }
        if(chlname_vec(i) == "LUT02"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(10);
        }
        if(chlname_vec(i) == "LUT03"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(11);
        }
        if(chlname_vec(i) == "XAT"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(12);
        }
        if(chlname_vec(i) == "XAT01"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(13);
        }
        if(chlname_vec(i) == "XAT02"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(14);
        }
        if(chlname_vec(i) == "XAT03"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(15);
        }
        if(chlname_vec(i) == "NEX"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(16);
        }
        if(chlname_vec(i) == "NEX01"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(17);
        }
        if(chlname_vec(i) == "NEX02"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(18);
        }
        if(chlname_vec(i) == "NEX03"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(19);
        }
        if(chlname_vec(i) == "CLAN"){
            dipolemat(0)(i) = dipolemat(0)(i) * -debye_vec(20);     //For HF the Minus goes Here
        }
        if(chlname_vec(i) == "CHLN"){
            dipolemat(0)(i) = dipolemat(0)(i) * debye_vec(21);     //For CAMB3LYP Here
        }
    }

    //cout<<dipolemat(0)<<endl;
    //cout<<"dipolemat:"<<dipolemat<<endl;   //errorchecking
    //cout<<"dipolemat:"<<dipolemat(0)<<endl;   //errorchecking
    //cout<<chlname_vec<<endl;
}

void get_nanc_nbnd(Vektor< Vektor< Vektor<double> > >& dipolemat){

    //At first let us obtain the unitvectors of the dipole transitionmoments from the NB-ND file
    //they will be stored as vectors within an Vektor

    Vektor< Vektor<double> > ei_dipole(0);

    ifstream nbndfile ("na_nc_nb_nd.dat");
    string data_str;
    string tmpstr;

    int k=0;

    Vektor<double> NA_dipole(3);
    Vektor<double> NB_dipole(3);
    Vektor<double> NC_dipole(3);
    Vektor<double> ND_dipole(3);

    Vektor<double> NBND(3);
    Vektor<double> NCNA(3);

    if (nbndfile.is_open()){
        while (nbndfile.good()){
        	if(nbndfile){

        	nbndfile >> tmpstr >> tmpstr >> tmpstr >> NA_dipole(0) >> NA_dipole(1) >> NA_dipole(2);
        	nbndfile >> tmpstr >> tmpstr >> tmpstr >> NB_dipole(0) >> NB_dipole(1) >> NB_dipole(2);
        	nbndfile >> tmpstr >> tmpstr >> tmpstr >> NC_dipole(0) >> NC_dipole(1) >> NC_dipole(2);
        	nbndfile >> tmpstr >> tmpstr >> tmpstr >> ND_dipole(0) >> ND_dipole(1) >> ND_dipole(2);


        	NBND = ND_dipole - NB_dipole;
        	NCNA = NA_dipole - NC_dipole;

        	ei_dipole.resize((int)ei_dipole.size() + 2);

        	ei_dipole(k) 		= NBND / norm_2(NBND);
        	ei_dipole(k + 1)	= NCNA / norm_2(NCNA);

            k = k + 2;
            }
            else break;
        }
    }

    //-----Kill Last 2 Elements (douplicates)-----
    ei_dipole.resize((int)ei_dipole.size() - 2);


    //-----Output the ei_dipoles-----
//	for(int i = 0; i < (int)ei_dipole.size(); i++){
//
//		cout <<"i = " << i << "  " << ei_dipole(i) << endl;
//	}

    //-----Dipolematrix: (0) UnitTransition NBND-vectors (1) Chl-Centers (2) NBND / NCNA Vectors----
    dipolemat.resize((int)dipolemat.size() + 1);

    dipolemat(2) = ei_dipole;



	//-----Method without clearing last 2 Elements
//				getline (nbndfile,data_str);
//
//	        	if(nbndfile){
//
//	        	stringstream data_stream(data_str);
//
//	        	data_stream >> tmpstr >> tmpstr >> tmpstr >> NA_dipole(0) >> NA_dipole(1) >> NA_dipole(2);
//
//	        	getline (nbndfile,data_str);
//	        	data_stream.clear();
//	        	data_stream << data_str;
//	        	data_stream >> tmpstr >> tmpstr >> tmpstr >> NB_dipole(0) >> NB_dipole(1) >> NB_dipole(2);
//
//	        	getline (nbndfile,data_str);
//	        	data_stream.clear();
//	        	data_stream << data_str;
//	        	nbndfile >> tmpstr >> tmpstr >> tmpstr >> NC_dipole(0) >> NC_dipole(1) >> NC_dipole(2);
//
//	        	getline (nbndfile,data_str);
//	        	data_stream.clear();
//	        	data_stream << data_str;
//	        	data_stream >> tmpstr >> tmpstr >> tmpstr >> ND_dipole(0) >> ND_dipole(1) >> ND_dipole(2);



    //cout<<ei_dipole<<endl;          //errorchecking
    //cout<<RCenter_dipole<<endl;     //errorchecking

    //Dipolematrix has UnitTransition Dipolevectors in & The Centers of the CHL for RM
//    dipolemat.resize(dipolemat.size()+2);
//
//    dipolemat(0)=ei_dipole;
//    dipolemat(1)=RCenter_dipole;

    //sets the dipolemoments according to init.dat debye_a, debye_b

    //cout<<"Dipolemat size: "<<dipolemat(0).size()<<endl;      //errorchecking
    //cout<<"Debye_Vec size: "<<debye_vec.size()<<endl;
    //cout<<"Chl_Namevec size: "<<chlname_vec.size()<<endl;
    //cout<<chlname_vec<<endl;
    //cout<<dipolemat(0)<<endl;
}

void turn_dipoles(Vektor< Vektor< Vektor<double> > >& dipolemat, double t_ang){

	Vektor<double> Turn_vec;
	double Pi = 3.14159265;

	for(int i = 0; i < (int)dipolemat(0).size(); i++){

		Turn_vec = cos( t_ang * Pi / 180.0 ) * dipolemat(2)(i * 2) + sin( t_ang * Pi / 180.0) *  dipolemat(2)(i * 2 + 1);
		Turn_vec = Turn_vec / norm_2(Turn_vec);

		Turn_vec = Turn_vec * norm_2( dipolemat(0)(i) );

		dipolemat(0)(i) = Turn_vec;
	}

}

void get_coordinates(Vektor< ChlAtoms >& Chl_Vec, string FileName){

	//-----Save the Coordinates of Chls into a Class-Vector-----

	ifstream xyz_file(FileName);
    string str;
    string tmp_string;

    Vektor< double > tmp_vec(3);

    int Number_Mols = 0;
    int Number_Coordinates = 0;

    //-----Get the # of Molecules and # of Coordinates in the Combined .xyz file-----
    if (xyz_file.is_open()) {

    		getline(xyz_file, str);
    		stringstream(str) >> Number_Mols;
    }

    Chl_Vec.resize(Number_Mols);

    //-----Save the Coordinates into the Class-Vector-----
    for(int k = 0; k < Number_Mols; k++){

    	getline(xyz_file, str);
    	stringstream(str) >> Number_Coordinates;

    	//-----Comment in xyz-Files-----
    	getline(xyz_file, str);

		for(int i = 0; i < Number_Coordinates; i++){

			getline(xyz_file, str);
			stringstream Coords(str);

			Coords >> tmp_string;

			for(int j = 0; j < 3; j++){
				Coords >> tmp_vec(j);
			}

			Chl_Vec(k).setAtomcoords(tmp_vec);
		}
    }
}

int get_charges(ChlAtoms& Chl_Vec, string FileName){

	ifstream xyz_file(FileName);
    string str;
    string tmp_string;
    double tmp_charge = 0;

    //-----Set the End of the File to the Line beginning with "Total"-----
    string stop_string = "Total";

    //-----Skip the first 2 lines in the fitcharges.out file-----
    if (xyz_file.is_open()) {

    		getline(xyz_file, str);
    		getline(xyz_file, str);
    }

    if (xyz_file.is_open()) {
    	while (xyz_file.good()) {
    		getline(xyz_file, str);

    		if (xyz_file) {

    			stringstream ChargeStream(str);

    			for (int i = 0; i < 5; i++) {

    				ChargeStream >> tmp_string;

    				if(tmp_string == stop_string){
//    					cout << "Breaking!" << endl;
    					return 0;
    				}
    			}

    			ChargeStream >> tmp_charge;

    			Chl_Vec.addChargeQy(tmp_charge);

    			}
    			else
    				break;
    	}
    }

    return 0;
}

int get_BChl_charges(ChlAtoms& Chl_Vec, string FileName){

	ifstream xyz_file(FileName);
    string str;
    string tmp_string;
    double tmp_charge = 0;

    //-----Real till EOF-----

	if (xyz_file.is_open()) {
		while (xyz_file.good()) {
			getline(xyz_file, str);

			if (xyz_file) {

				stringstream ChargeStream(str);

//    			for (int i = 0; i < 5; i++) {
//
//    				ChargeStream >> tmp_string;
//
//    				if(tmp_string == stop_string){
////    					cout << "Breaking!" << endl;
//    					return 0;
//    				}
//    			}

				ChargeStream >> tmp_charge;

				Chl_Vec.addChargeQy(tmp_charge);

			}
			else
				break;
		}
	}

    return 0;
}

void get_coordinates_pdb(Vektor< ChlAtoms >& Chl_Vec, string FileName){

    ifstream PDB;
    PDB.open(FileName);

    //-----EINLESEN DER DATEN ANHAND DES PDB-FILES-----

    string ATOM;
	string AtomType;
	string resType;
	int resNumber = 0;

	Vektor<double> coord(3);

    string str;
    string tmp_str;

    if (PDB.is_open()) {
		while (PDB.good()) {
			getline(PDB, str);
			if (PDB) {
				stringstream sstest(str);

				sstest >> ATOM;

				if (ATOM == "HETATM") {

					sstest >> tmp_str;
					sstest >> AtomType;
					sstest >> resType;
					sstest >> resNumber;

					if (AtomType == "MG" && (resType == "CLA" || resType == "CHL")) {

						Chl_Vec.resize(Chl_Vec.size() + 1);
						Chl_Vec(Chl_Vec.size() - 1).setAtomtype(resType);
						Chl_Vec(Chl_Vec.size() - 1).resNumber = resNumber;
					}

					Chl_Vec(Chl_Vec.size() - 1).setAtomname(AtomType);

					//-----Get to the Coordinates-----
					for (int i = 0; i < 3; i++) {
						sstest >> coord(i);
					}

					Chl_Vec(Chl_Vec.size() - 1).setAtomcoords(coord);
				}
			}
		}
    }

    PDB.close();
}

void setup_LHCII(Vektor< ChlAtoms >& Chl_Vec){

	//-----Setup the LHCII Trimer-----

    get_coordinates_pdb(Chl_Vec, "trimer_chl_rotated.pdb");

    //-----Löschen des Phetylrestes aus dem Pigments Klassen-Vektor-----
    //-----Löscht die Einträge C2, C3, usw. aus den Coordinaten und Atomnamen der Chl-----
	Vektor<string> phetylrest(19);
	for (int i = 0; i < (int) phetylrest.size(); i++) {
		phetylrest(i) = "C" + to_string(i + 2);
	}

	for (int i = 0; i < (int) Chl_Vec.size(); i++) {

		for (int j = 0; j < (int) phetylrest.size(); j++) {
			Chl_Vec(i).deletefromlist(phetylrest(j));
		}
	}

	//-----Read in the Chla/b Charges - The SAME routine is used as for BChl bc of the SAME File-Format!-----
    for(int h = 0; h < (int)Chl_Vec.size(); h++){

    	if(Chl_Vec(h).atomtype == "CLA"){
    		get_BChl_charges(Chl_Vec(h), "chla_charges_B3LYP.dat");
    	}

    	else if(Chl_Vec(h).atomtype == "CHL"){
    		get_BChl_charges(Chl_Vec(h), "chlb_charges_B3LYP.dat");
    	}

    	else
    		cout << "ERROR: Charge-File not Found!" << endl;
    }
}

void setup_LHCII_monomer(Vektor< ChlAtoms >& Chl_Vec){

	//-----Setup the LHCII Trimer-----

    get_coordinates_pdb(Chl_Vec, "monomer_chl.pdb");

    //-----Löschen des Phetylrestes aus dem Pigments Klassen-Vektor-----
    //-----Löscht die Einträge C2, C3, usw. aus den Coordinaten und Atomnamen der Chl-----
	Vektor<string> phetylrest(19);
	for (int i = 0; i < (int) phetylrest.size(); i++) {
		phetylrest(i) = "C" + to_string(i + 2);
	}

	for (int i = 0; i < (int) Chl_Vec.size(); i++) {

		for (int j = 0; j < (int) phetylrest.size(); j++) {
			Chl_Vec(i).deletefromlist(phetylrest(j));
		}
	}

	//-----Read in the Chla/b Charges - The SAME routine is used as for BChl bc of the SAME File-Format!-----
    for(int h = 0; h < (int)Chl_Vec.size(); h++){

    	if(Chl_Vec(h).atomtype == "CLA"){
    		get_BChl_charges(Chl_Vec(h), "chla_charges_B3LYP.dat");
    	}

    	else if(Chl_Vec(h).atomtype == "CHL"){
    		get_BChl_charges(Chl_Vec(h), "chlb_charges_B3LYP.dat");
    	}

    	else
    		cout << "ERROR: Charge-File not Found!" << endl;
    }
}

void setup_ChlVec(Vektor< ChlAtoms >& Chl_Vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	//-----Calculate Chl Properties for every Chl that is Read-In-----

	for (int i = 0; i < (int) Chl_Vec.size(); i++) {

//		Chl_Vec(i).calculate_Dipole_Moment();

		Chl_Vec(i).set_Dipole_Moment_NBND();

		if(Chl_Vec(i).atomtype == "CLA"){
			Chl_Vec(i).rescale_Dipole_Moment(4.0);
		}
		if(Chl_Vec(i).atomtype == "CHL"){
			Chl_Vec(i).rescale_Dipole_Moment(3.4);
		}

		Chl_Vec(i).set_Mol_Center();

		//cout << "MolCenter from Positions: " << Chl_Vec(i).atomcoords(5) << "   " << Chl_Vec(i).atomcoords(33) << endl << endl;

		Chl_Vec(i).calculate_Quadrupole_Moment();
		Chl_Vec(i).calculate_Intrinsic_Quadrupole_Moment();

//		cout << "Dipole Moment: " << Chl_Vec(i).Dipole_Mom << " Center: " << Chl_Vec(i).Mol_Center << endl;
//		cout << "Quadrupole Moment: " << Chl_Vec(i).Quadrupole_Moment << endl;
//		cout << "Intrinsic Quadrupole Moment: " << Chl_Vec(i).Intrinsic_Quadrupole_Moment << endl;
//

		//-----This is done for the old Structure of the Program-----
		dipolemat(0)(i) = Chl_Vec(i).Dipole_Mom;
		dipolemat(1)(i) = Chl_Vec(i).Mol_Center;

	}
}

void GetCoordinates_xyz(Vektor< ChlAtoms >& Chl_Vec, string FileName, long numb_mol){

	//-----Save the Coordinates of Chls into a Class-Vector-----

	ifstream xyz_file(FileName);
    string str;
    string tmp_string;

    Vektor< double > tmp_vec(3);

    int numb_coords = 0;

    //-----Get the # of Molecules and # of Coordinates in the Combined .xyz file-----
    if (xyz_file.is_open()) {

    		getline(xyz_file, str);
    		stringstream(str) >> numb_coords;
    }
    else{
    	cout << "File '" << FileName << "' not found in Directory!" << endl;
    }
    //-----Comment in xyz-Files-----
    getline(xyz_file, str);
    Chl_Vec.resize(numb_mol);

    //-----Save the Coordinates into the Class-Vector-----
    for(int k = 0; k < numb_mol; k++){

		for(int i = 0; i < numb_coords / numb_mol; i++){

			getline(xyz_file, str);
			stringstream Coords(str);

			Coords >> tmp_string;

			for(int j = 0; j < 3; j++){
				Coords >> tmp_vec(j);
			}

			Chl_Vec(k).setAtomname(tmp_string);
			Chl_Vec(k).setAtomcoords(tmp_vec);
		}
    }
}

void setup_Helix(Vektor< ChlAtoms >& Chl_Vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	long numb_mol = 16;

//	GetCoordinates_xyz(Chl_Vec, "SQ2_Structure_reordered.xyz", numb_mol);
//	GetCoordinates_xyz(Chl_Vec, "SQ34_Check.xyz", numb_mol);
//	GetCoordinates_xyz(Chl_Vec, "SQ1-3_Structure.xyz", numb_mol);
//	GetCoordinates_xyz(Chl_Vec, "SQ8_Structure_reordered.xyz", numb_mol);
//	GetCoordinates_xyz(Chl_Vec, "SQ16_Structure.xyz", numb_mol);

	//-----Check if Helix is Left-Handed-----
//	GetCoordinates_xyz(Chl_Vec, "SQ16_left_handed_helix.xyz", numb_mol);
	GetCoordinates_xyz(Chl_Vec, "SQ16_xyz_lh_squeezed_helix.dat", numb_mol);

	int i = 0;

	for(auto& Chl : Chl_Vec){
		get_charges(Chl, "fitted_charges_Trans.out");
		Chl.calculate_Dipole_Moment();
//		Chl.rescale_Dipole_Moment(10.0);
//		Chl.set_Mol_Center(11, 13);
		Chl.set_Mol_Center(14, 14);

//		if(i == 1 || i == 4 || i == 6){
//			Chl.Dipole_Mom = -Chl.Dipole_Mom;
//		}

		cout << Chl.Dipole_Mom << endl;

		//-----This is done for the old Structure of the Program-----
		dipolemat(0)(i) = Chl_Vec(i).Dipole_Mom;
		dipolemat(1)(i) = Chl_Vec(i).Mol_Center;
		i++;
	}
}

void setup_Helix_2StateModel(Vektor< ChlAtoms >& Chl_Vec, Vektor< Vektor< Vektor<double> > >& dipolemat){

	long numb_mol = 16;

	GetCoordinates_xyz(Chl_Vec, "SQ16_xyz_lh_helix.dat", numb_mol);

	int i = 0;
	for(auto& Chl : Chl_Vec){
		get_charges(Chl, "fitted_charges_Trans.out");
		Chl.calculate_Dipole_Moment();
//		Chl.rescale_Dipole_Moment(10.0);
//		Chl.set_Mol_Center(11, 13);
		Chl.set_Mol_Center(14, 14);

		cout << Chl.Dipole_Mom << endl;

		//-----This is done for the old Structure of the Program-----
		dipolemat(0)(i) = Chl_Vec(i).Dipole_Mom;
		dipolemat(1)(i) = Chl_Vec(i).Mol_Center;
		i++;
	}

	Vektor< ChlAtoms > Chl_Vec_XS2;
	GetCoordinates_xyz(Chl_Vec_XS2, "SQ16_xyz_lh_helix.dat", numb_mol);

	for(auto& Chl : Chl_Vec_XS2){
		get_charges(Chl, "fitted_charges_Trans_XS2.out");
		Chl.calculate_Dipole_Moment();

		Chl.set_Mol_Center(14, 14);

		cout << Chl.Dipole_Mom << endl;

		//-----This is done for the old Structure of the Program-----
		dipolemat(0)(i) = Chl_Vec_XS2(i - numb_mol).Dipole_Mom;
		dipolemat(1)(i) = Chl_Vec_XS2(i - numb_mol).Mol_Center;
		i++;
	}

	for(auto Chl : Chl_Vec_XS2){
		push_back(Chl_Vec, Chl);
	}
}

#endif // INIT_HPP_INCLUDED
