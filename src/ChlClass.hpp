#ifndef CHLCLASS_HPP_INCLUDED
#define CHLCLASS_HPP_INCLUDED
#include <string>
#include "Utilities.hpp"
using namespace std;


class ChlAtoms{

public:

    Vektor<Vektor<double> > atomcoords;
    Vektor<Vektor<double> > EDatomcoordsBY;
    Vektor<Vektor<double> > EDatomcoordsBX;
    Vektor<string> atomname;
    Vektor<double> chargesQy;
    Vektor<double> chargesBy;
    Vektor<double> chargesBx;
    string atomtype;
    int resNumber;

    Vektor< double > Dipole_Mom;
    Matrix< double > Quadrupole_Moment;
    Vektor< double > Mol_Center;
    Matrix< double > Intrinsic_Quadrupole_Moment;

    //Vektor<Atoms> atom;


    void setAtomname(string AtName){
        atomname.resize(atomname.size()+1);
        atomname(atomname.size()-1) = AtName;
    }

    //*****Set Charges for Qy, By Bx; Set Atomcoords or whole Atomcoords Vector, Set Atomtype or delete cetrain amount of Atoms*****

    void addChargeQy(double inputchargesQy){
    	chargesQy.resize((int)chargesQy.size() + 1);
    	chargesQy((int)chargesQy.size() - 1) = inputchargesQy;
    }

    void setChargesQy(Vektor<double> inputchargesQy){

        chargesQy = inputchargesQy;
    }

    void setChargesBy(Vektor<double> inputchargesBy){

        chargesBy = inputchargesBy;
    }

    void setChargesBx(Vektor<double> inputchargesBx){

        chargesBx = inputchargesBx;
    }

    void setAtomcoords(Vektor<double> coordsofoneatom){

        atomcoords.resize(atomcoords.size()+1);
        atomcoords(atomcoords.size()-1) = coordsofoneatom;
    }

    void setAtomcoordsVector(Vektor<Vektor<double> > inputcoords){

        atomcoords = inputcoords;
    }

    void setEDAtomcoordsBY(Vektor<Vektor<double> > inputcoords){

        EDatomcoordsBY = inputcoords;
    }

    void setEDAtomcoordsBX(Vektor<Vektor<double> > inputcoords){

        EDatomcoordsBX = inputcoords;
    }

    void setAtomtype(string type){

        atomtype = type;
    }

    void deletefromlist(string AtName){


        int isentry = 0;
		int position = 0;

		for (int i = 0; i < (int) atomname.size(); i++) {

			if (atomname(i) == AtName) {
				isentry = 1;
				position = i;
			}
		}

		if (isentry == 1) {
			for (int j = position; j < (int) atomname.size() - 1; j++) {

				atomname(j) = atomname(j + 1);
				atomcoords(j) = atomcoords(j + 1);

			}

			atomname.resize((int) atomname.size() - 1);
			atomcoords.resize((int) atomcoords.size() - 1);
		}

		else if (isentry == 0) {
//			cout << "ERROR: void deletefromlist() element was not found. This is Chl " << resNumber << " reporting!" << endl;
		}

    }

    void calculate_Dipole_Moment(){

    	Dipole_Mom.resize(3);
    	for(auto& i : Dipole_Mom){
    		i = 0;
    	}

    	for(int i = 0; i < (int)atomcoords.size(); i++){

    		Dipole_Mom = Dipole_Mom + atomcoords(i) * chargesQy(i) * 4.8;
    	}
    }

    void set_Dipole_Moment_NBND(){

    	Vektor<double> NB_vec(3);
    	Vektor<double> ND_vec(3);

    	for(int i = 0; i < (int)atomname.size(); i++){
    		if(atomname(i) == "NB"){
    			NB_vec = atomcoords(i);
    		}
    		if(atomname(i) == "ND"){
    			ND_vec = atomcoords(i);
    		}
    	}

		Dipole_Mom.resize(3);
		Dipole_Mom = ND_vec - NB_vec;
	}

    void rescale_Dipole_Moment(double Dip_Strength){

    	double Calc_Dip_Strength = 0;
    	Calc_Dip_Strength = norm_2(Dipole_Mom);

//    	chargesQy = chargesQy / Calc_Dip_Strength * Dip_Strength;
    	Dipole_Mom = Dipole_Mom / Calc_Dip_Strength * Dip_Strength;
    }

    void set_Mol_Center(Vektor<double> Center){
    	Mol_Center = Center;
    }

    void set_Mol_Center(double pos1, double pos2){
        	Mol_Center = (atomcoords(pos1) + atomcoords(pos2)) / 2.0;
    }

    void set_Mol_Center(){

    	Vektor<double> NB_vec(3);
    	Vektor<double> ND_vec(3);

    	for(int i = 0; i < (int)atomname.size(); i++){
    		if(atomname(i) == "NB"){
    			NB_vec = atomcoords(i);
    		}
    		if(atomname(i) == "ND"){
    			ND_vec = atomcoords(i);
    		}
    	}

    	Mol_Center.resize(3);
    	Mol_Center = (ND_vec + NB_vec) / 2.0;
	}

    void calculate_Quadrupole_Moment(){

//    	Vektor<Vektor<double> > relative_atomcoords;
//    	relative_atomcoords = atomcoords;
//
//    	for(int i = 0; i < (int)atomcoords.size(); i++){
//    		relative_atomcoords(i) = atomcoords(i) - Mol_Center;
//    	}

    	Quadrupole_Moment.resize(3,3);

    	for(int i = 0; i < 3; i++){
    		for(int j = 0; j < 3; j++){

    			Quadrupole_Moment(i,j) = 0.0;
    		}
    	}


    	for(int i = 0; i < 3; i++){
    		for(int j = 0; j < 3; j++){

    			for(int k = 0; k < (int)atomcoords.size(); k++){

    				Quadrupole_Moment(i,j) = Quadrupole_Moment(i,j) + atomcoords(k)(i) * atomcoords(k)(j) * chargesQy(k) * 4.8;
    			}

//    			cout << " Quadrupole Moment i j : " << i << " " << j << "  " << Quadrupole_Moment << endl;
    		}
    	}

    }

    void calculate_Intrinsic_Quadrupole_Moment(){

        	Vektor<Vektor<double> > relative_atomcoords;
        	relative_atomcoords = atomcoords;

        	for(int i = 0; i < (int)atomcoords.size(); i++){
        		relative_atomcoords(i) = atomcoords(i) - Mol_Center;
        	}

        	Intrinsic_Quadrupole_Moment.resize(3,3);

        	for(int i = 0; i < 3; i++){
        		for(int j = 0; j < 3; j++){

        			Intrinsic_Quadrupole_Moment(i,j) = 0.0;
        		}
        	}


        	for(int i = 0; i < 3; i++){
        		for(int j = 0; j < 3; j++){

        			for(int k = 0; k < (int)atomcoords.size(); k++){

        				Intrinsic_Quadrupole_Moment(i,j) = Intrinsic_Quadrupole_Moment(i,j) + relative_atomcoords(k)(i) * relative_atomcoords(k)(j) * chargesQy(k) * 4.8;
        			}

    //    			cout << " Quadrupole Moment i j : " << i << " " << j << "  " << Quadrupole_Moment << endl;
        		}
        	}

        }

    void Rotate_Dipole_Moment(double ang) {

    	Vektor<double> NA_vec(3);
		Vektor<double> NB_vec(3);
		Vektor<double> NC_vec(3);
		Vektor<double> ND_vec(3);

		for (int i = 0; i < (int) atomname.size(); i++) {

			if (atomname(i) == "NA") {
				NA_vec = atomcoords(i);
			}
			if (atomname(i) == "NB") {
				NB_vec = atomcoords(i);
			}
			if (atomname(i) == "NC") {
				NC_vec = atomcoords(i);
			}
			if (atomname(i) == "ND") {
				ND_vec = atomcoords(i);
			}
		}

		Vektor<double> Turn_vec(3);
		Vektor<double> NBND_vec(3);
		Vektor<double> NANC_vec(3);

		NBND_vec = ND_vec - NB_vec;
		NANC_vec = NC_vec - NA_vec;

		double Pi = 3.14159265;

		Turn_vec = cos( ang * Pi / 180.0 ) * NANC_vec + sin( ang * Pi / 180.0) * NBND_vec;
		Turn_vec = Turn_vec / norm_2( Turn_vec ) * norm_2( Dipole_Mom );
//		Turn_vec = Turn_vec * norm_2( Dipole_Mom );
//		Dipole_Mom.resize(3);

		Dipole_Mom = Turn_vec;
	}


};

#endif // CHLCLASS_HPP_INCLUDED
