/*
 *  Created by Ryky Nelson Aug 2019
 */

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <sstream>
#include <iomanip>
#include <boost/filesystem.hpp>
#include <boost/format.hpp>
#include <cstring>
#include <cstdlib>
#include "elements.h"

using namespace std;
using namespace Eigen;
using namespace boost::filesystem;

void Usage(string exe_name){    
    cerr << "please specify your POSCAR or CONTCAR;\n"
	 << "e.g.: " << exe_name << " POSCAR\n";
}

void author(){    
    cout << "Program is created by Ryky Nelson.\n";
}

string GetEnv( const string & var ) {
     const char * val = ::getenv( var.c_str() );
     if ( val == 0 ) {
         return "";
     }
     else {
         return val;
     }
}

int main(int argc, char* argv[]){

    if (argc < 2 ) {
	Usage(argv[0]);
	return 1;
    }

    string flag(argv[1]);
    if (flag == "-c" || flag == "-C") {
	author();
	return 0;
    }

    string vaspin = "";
    vaspin = argv[1];

    std::ifstream vaspFile( vaspin.c_str() );
    if( !vaspFile || (vaspFile.peek() == std::ifstream::traits_type::eof()) ) {
    	cerr << "Unable to open " << vaspin << " or it is empty!" << endl;
    	return 0;
    }
    
    string line; stringstream ss;
    string trash, sval;
    double dval;

    getline(vaspFile,line); // ignore the header

    getline(vaspFile,line);
    ss.str(line); ss.clear();
    double acell;
    ss >> acell;

    // get the primitive vector from POSCAR 
    Matrix3d primVec = Matrix3d::Zero();
    for (int idummy = 0; idummy < 3; ++idummy) { // the elements' order is kept!!!
	getline(vaspFile,line);
	ss.str(string()); ss.clear(); ss.str(line);
	ss >> primVec(idummy,0) >> primVec(idummy,1) >> primVec(idummy,2);
    }
    double maxEl = primVec.array().abs().maxCoeff();
    acell *= maxEl;
    primVec /= maxEl;

    // get the atoms and position-vectors from POSCAR
    size_t natom(0), ntyp(0);
    string atyp;
    vector<string> atomType;
    getline(vaspFile,line);
    ss.str(string()); ss.clear(); ss.str(line);
    while(ss>>atyp) atomType.push_back(atyp);
    ntyp = atomType.size();

    vector<size_t> numPerType;
    getline(vaspFile,line);
    ss.str(string()); ss.clear(); ss.str(line);
    for(vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
	size_t tdummy; ss >> tdummy;
	natom += tdummy;
	numPerType.push_back(tdummy);
    }
    
    getline(vaspFile,line); // kind of coordinates: D or C, for now assuming D

    vector<Vector3d> apos;
    for(size_t idummy = 0; idummy < ntyp; ++idummy) {
	size_t npTyp = numPerType[idummy];
	for(size_t jdummy = 0; jdummy < npTyp; ++jdummy) {
	    getline(vaspFile,line);
	    ss.str(string()); ss.clear(); ss.str(line);
	    Vector3d pos;
	    ss >> pos(0) >> pos(1) >> pos(2);
	    apos.push_back(pos);
	}
    }

    // print output
    cout << "acell " << "3*" << acell << " angstrom\n";
    // print rprim
    cout << fixed << setprecision(12);
    cout << "\nrprim\n" << primVec << endl
	 << "\nnatom " << natom << "\n\nntypat " << ntyp << endl
	 << "\ntypat "; 
    for(size_t idummy = 0; idummy < ntyp; ++idummy)
	for(size_t jdummy = 0; jdummy < numPerType[idummy]; ++jdummy) cout << (idummy+1) << ' ';

    cout << "\n\nznucl ";
    for (vector<string>::iterator it = atomType.begin(); it != atomType.end(); ++it) {
	int atomOfE = elementsMass.at(*it);
	cout << atomOfE << ' ';
    }
    cout << endl << endl;

    // print xred
    cout << "xred\n" << fixed << setprecision(12);
    size_t iterdummy(0);
    for(size_t idummy = 0; idummy < ntyp; ++idummy) {
	size_t npTyp = numPerType[idummy];
	for(size_t jdummy = 0; jdummy < npTyp; ++jdummy) {
	    cout << ' ' << apos[iterdummy].transpose() << endl;
	    ++iterdummy;
	}
    }

    cout << "\necut XXECUTXX\n"
	 << "\npawecutdg XXPAWECUTDGXX\n"
	 << "\nnband XXNBNDXX\n";

    cout << "\ntolvrs 1.0d-12\n"
	 << "\nnstep 100\n"
	 << "\nistwfk *1\n";

    // string nspin("1");
    string NSpin = GetEnv("NSPIN");
    if ( NSpin.empty() ) NSpin = '1';
    ss.str(string()); ss.clear(); ss.str(NSpin);
    int nspin;
    ss >> nspin;
    cout << "\nnsppol " << nspin << endl;
    if (nspin == 2) {
	string magmom = GetEnv("MAGMOM");
	ss.str(string()); ss.clear(); ss.str(magmom);

	cout << "\nspinat\n" << fixed << setprecision(3);
	size_t tdummy;
	while(ss>>tdummy) {
	    double magmomV; ss >> magmomV;
	    for ( size_t idummy = 0; idummy < tdummy; ++idummy ) {
		cout << " 0.000 0.000 " << magmomV << endl;
	    }
	}
    }
	
    cout << "\noccopt 6\n";

    // print K_POINTS {automatic}
    cout << "\nngkpt ";
    string nkptString;
    int kptDensity(-1);
    istringstream itmp;
    if (argc > 2) itmp.str(argv[2]);
    if ( !(itmp >> kptDensity) || ! itmp.eof()) kptDensity = -1;       
    if (kptDensity > 0){
	int ngrid = kptDensity / natom;
	int mult = pow(( ngrid 
			* primVec.row(0).norm()
			* primVec.row(1).norm() 
			* primVec.row(2).norm() ), 1.0/3.0);
	
	Array3i nkpt;
	for (int ikpt = 0; ikpt < 3; ++ikpt) {
	    nkpt(ikpt) = round(mult / primVec.row(ikpt).norm());
	    if (nkpt(ikpt) == 0) nkpt(ikpt) = 1;
	}
	cout << nkpt.transpose() << endl;
    }
    else cout << "XXKPOINTXX\n";
    cout << "\nshiftk 0.0 0.0 0.0\n";


    return 0;

}
