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
  Vector3d rotvT; // inverted rotational vector
  Vector3d ang; // grain rotation Euler angles
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
    bool          init;
    MatrixXd      cell;

    CONFIG(Vector3d l, double ac, unsigned grains, unsigned total, MatrixXd cell) : l(l), ac(ac), grains(grains), grain(grains), atom_box(0), cell(cell)
    {
    atoms_box = 0;
    atoms_grain = pow(ceil(pow(pow(3.0,3.0/2.0) * total / pow(grains,1.0/3.0) / cell.cols(), 1.0/3.0)),3.0) * cell.cols();
    atom_grain.resize(atoms_grain);
    atom_box.resize(2*atoms_grain);
    shift = l/2.0;
    init = false;
    }
    
};  

#endif
