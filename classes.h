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
    int           grains;     // number of grains
    vector<GRAIN> grain;      // coordinates and properties of grain
    vector<ATOM>  atom_box;   // coordinates of atoms in the box
    int           atoms_box;  // total number of atoms in the box
    vector<ATOM>  atom_grain; // atoms in a current grain
    int           atoms_grain;// atoms in grain
    bool          init;       // init random generator
    MatrixXd      cell;       // unit cell vectors
    int           cell_type;  // unit cell code
    int           out_type;   // output file code
    double        v;          // volume per atom
    string        name;       // element name

    CONFIG(Vector3d l, double ac, unsigned grains, double v, MatrixXd cell) : l(l), ac(ac), grains(grains), grain(grains), atom_box(0), cell(cell), v(v)
    {
    atoms_box = 0;
    atoms_grain = pow(ceil(pow(pow(3.0,3.0/2.0) * ceil(l.prod()/v) / pow(grains,1.0/3.0) / cell.cols(), 1.0/3.0)),3.0) * cell.cols();
    atom_grain.resize(atoms_grain);
    atom_box.resize(2*atoms_grain);
    shift = l/2.0;
    init = false;
    }
    
};  

#endif
