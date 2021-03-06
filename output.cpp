#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <time.h>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#ifndef  CLASSES_H  
#include "classes.h"
#endif

using namespace Eigen;
using namespace std;

unsigned WriteLAMMPS(CONFIG &config)
{
    cout << "Saving LAMMPS config with " << config.atoms_box << " atoms." << endl;
    stringstream filename;
    string filename_s;


    filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains << ".in";
    filename_s = filename.str();
    ofstream outfile (filename_s.c_str(), ios::trunc);
    outfile.precision(5);

    Vector3d l;

    outfile << "Polycrystal made by polymaker. Box " << config.l.transpose() << " A, " << config.v << " A^3 per atom." << endl;
    outfile << "1 atom types\n";
    outfile << config.atoms_box << " atoms\n\n";
    outfile << -config.shift(0) << " " << config.shift(0) << " xlo xhi\n";
    outfile << -config.shift(1) << " " << config.shift(1) << " ylo yhi\n";
    outfile << -config.shift(2) << " " << config.shift(2) << " zlo zhi\n\n";
    //outfile << "Pair Coeffs\n\n";
    //outfile << "2 0.42276 4.099\n";
    outfile << endl;
    outfile << "Masses\n\n";
    outfile << "1 XXX\n";
    //outfile << "2 60\n";
    //outfile << "3 131.2937\n";
    outfile << endl;
    outfile << "Atoms\n\n";

    //pair_style      buck 8.0
    //pair_coeff      1 1 177079.254302103 0.343846  6917.93499044

    for (int i = 0; i < config.atoms_box ; i++)
    {
        outfile << i+1 << " 1 " << setw (15)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << endl;
    }
    
    outfile.close();

    return 0;
}

unsigned WriteXYZ(CONFIG &config)
{
    stringstream filename;
    string filename_s;

    filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains   << ".xyz";
    filename_s = filename.str();
    ofstream outfile (filename_s.c_str(), ios::trunc);
    outfile.precision(5);

    outfile << config.atoms_box << endl << "Polycrystal made by polymaker. Box " << config.l.transpose() << " A, " << config.v << " A^3 per atom." << endl;


    for (int  i = 0; i < config.atoms_box; i++)
    {
        outfile << config.name << setw (15)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << " " << config.atom_box[i].type << endl;
    }
    
    outfile.close();

    return 0;
}

unsigned WriteDLPOLY(CONFIG &config)
{
    stringstream filename;
    string filename_s;

    filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains   << "";
    filename_s = filename.str();
    ofstream outfile (filename_s.c_str(), ios::trunc);
    outfile.precision(5);

    outfile << "Polycrystal made by polymaker. Box " << config.l.transpose() << " A, " << config.v << " A^3 per atom." << endl;
    outfile << "0         3" << endl;
    outfile << config.l(0) << " 0.0000000000 " << " 0.0000000000"  <<endl;
    outfile << "0.0000000000 " << config.l(1) << " 0.0000000000"  <<endl;
    outfile << "0.0000000000 " << " 0.0000000000 "  << config.l(2) <<endl;


    for (int  i = 0; i < config.atoms_box; i++)
    {
        outfile << config.name << "   " << i+1 << endl;
        outfile << setw (22)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << endl;
    }
    
    outfile.close();

    return 0;
}

unsigned WriteGrains(CONFIG &config)
{

    stringstream filename;
    string filename_s;

    filename << "GRAINS";
    filename_s = filename.str();
    ofstream outfile (filename_s.c_str(), ios::trunc);
    outfile.precision(20);

    for (int  i = 0; i < config.grains; i++)
    {
        outfile << setw (15) << i << " " << (config.grain[i].r.array()/config.l.array()).transpose() << " " << config.grain[i].angle.transpose() << endl;
    }
    
    outfile.close();


    return 0;
}