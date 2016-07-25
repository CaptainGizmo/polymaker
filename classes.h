#ifndef   CLASSES_H
#define   CLASSES_H

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

const int AtomsPerCell = 100;
using namespace Eigen;
using namespace std;

class ATOM { // class for each atom's properties
    public:
        unsigned type;
        Vector3d r;
        //int      cell; // number of cell id in the grid
};

class GRAIN// : virtual ATOM
{
    public:
        Vector3d r; // grain center positions
        Vector3d rotvT; // inverted rotational vector
        Vector3d angle; // grain rotation Euler angles
};

class GRID
{
    public:
        //Vector3d r;  // physical grid position
        thrust::host_vector<int> id; //atom's ids
        Vector3i r;
};

class CONFIG //: virtual GRAIN, virtual ATOM
{ // full config of simulation
    public:
        Vector3d l;  // size of box in A
        Vector3d shift;
        double ac;      // lattice const
        int           grains;     // number of grains
        vector<GRAIN> grain;      // coordinates and properties of grain
        vector<ATOM>  atom_box;   // coordinates of atoms in the box
        int           atoms_box;  // total number of atoms in the box
        vector<ATOM>  atom_grain; // atoms in a current grain
        int           atoms_grain;// atoms in grain
        bool          init;       // init random generator
        bool          time;       // running time output
        bool          read_grains;// read grains from the file
        MatrixXd      unit_cell;  // unit cell vectors
        int           cell_type;  // unit cell code
        int           out_type;   // output file code
        double         v;          // volume per atom
        string        name;       // element name
        string        cell;       // cell type name
        string        filename;   // grains parameter file name
        vector<GRID>  grid;       // nearest neighbors greed
        int           grid_size;  // total number of elements in grid = a^3

    CONFIG(Vector3d l, double ac, unsigned grains, double v, MatrixXd unit_cell) : l(l), ac(ac), grains(grains), grain(grains), unit_cell(unit_cell), v(v)
    {
        atoms_box = 0;
        atoms_grain = pow(ceil(pow(pow(3.0,3.0/2.0) * ceil(l.prod()/v) / pow(grains,1.0/3.0) / unit_cell.cols(), 1.0/3.0)),3.0) * unit_cell.cols();
        atom_grain.resize(atoms_grain);
        atom_box.resize(2*atoms_grain);
        grid_size = ceil( pow( atoms_grain/AtomsPerCell, 1.0/3.0 ));
        grid.resize( pow( grid_size, 3.0) );
        shift = l/2.0;
        init = false;
        time = false;
        read_grains = false;
    }

};

#endif
