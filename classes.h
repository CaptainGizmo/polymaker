#ifndef   CLASSES_H
#define   CLASSES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

using namespace Eigen;
using namespace std;

class ATOM { // class for each atom's properties
public:
    unsigned type;
    Vector3d r;
};

class GRAIN// : virtual ATOM
{
public:
  Vector3d r; // grain center positions
};


class CONFIG //: virtual GRAIN, virtual ATOM
{ // full config of simulation
public:
    Vector3d l;  // size of box in A
    Vector3d shift;
    float ac;      // lattice const
    //float fnn; // first nearest neigbour
    int           grains;    // number of grains
    vector<GRAIN> grain;
    vector<ATOM>  atom_box;
    int           atoms_box;
    vector<ATOM>  atom_grain;
    int           atoms_grain;

    CONFIG(Vector3d l_, double ac_, unsigned grains_, unsigned total_) : l(l_), ac(ac_), grains(grains_), grain(grains_), atom_box(2*total_)
    {
    atoms_box=0;
    atoms_grain=4*pow(ceil(pow(total_/4.0,1.0/3.0)),3);
    //atoms_grain=2*pow(ceil(pow(total_/2.0,1.0/3.0)),3);
    atom_grain.resize(atoms_grain);
    shift = l/2.0;
    }
    
};  

#endif
