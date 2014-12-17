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
#include <iterator>

#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/LU>

#include <boost/program_options.hpp>

#ifndef  CLASSES_H  
#include "classes.h"
#endif

using namespace Eigen;
using namespace std;
namespace po = boost::program_options;

template<class T>
ostream& operator<<(ostream& os, const vector<T>& v)
{
    copy(v.begin(), v.end(), ostream_iterator<T>(cout, " ")); 
    return os;
}

CONFIG getinput(int argc, char **argv)
{

    int grains, cell_type, out_type;
    double ac, v;
    string cell, name, out;
    MatrixXd unit_cell;
    Vector3d l;

    try
    {
      // Declare the supported options.
      po::options_description desc("Allowed options");
      desc.add_options()
          ("help", "produce help message")

          //("optimization",   po::value<int>(&opt)->default_value(10),    "optimization level")
          //("include-path,I", po::value< vector<string> >()->composing(), "include path")

          ("V,v", po::value<double>(&v)->required(), "volume per atom")
          ("N,n", po::value<int>(&grains)->default_value(1), "number of grains")
          ("A,a", po::value<double>(&l[0])->required(), "box side a")
          ("B,b", po::value<double>(&l[1]), "box side b")
          ("C,c", po::value<double>(&l[2]), "box side c")
          ("type,t", po::value< vector<string> >(), "lattice type")
          ("out,o", po::value< vector<string> >(), "output file type")
          ("element,e", po::value< vector<string> >(), "element name")
;
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);    

      if (vm.count("help")) {
        cout << desc << "\n";
        exit(1);
      }

      if (!vm.count("B")) l[1] = l[0];
      if (!vm.count("C")) l[2] = l[0];
      
      if (vm.count("out"))
      {
          cout << "Output file type is " << vm["out"].as< vector<string> >() << "\n";
      }
    
      name = "Xe";
      cell = "fcc";
      out = "xyz";
    
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        exit(1);
    } 
    
    if (cell.compare("bcc") == 0)
    {
        cell_type = 1;
    }
    else 
    if (cell.compare("fcc") == 0)
    {
       cell_type = 2;
    }
    else
    {
       throw "Unknown unit cell!";
    }

    if (out.compare("dlpoly") == 0)
    {
        out_type = 1;
    }
    else 
    if (out.compare("lammps") == 0)
    {
       out_type = 2;
    }
    else 
    if (out.compare("xyz") == 0)
    {
       out_type = 3;
    }
    else
    {
       throw "Unknown output format!";
    }



    try{
      switch(cell_type)
      {
        case 1: 
            // bcc
            unit_cell.resize(3,2);
            unit_cell << 0, 0.5,
                        0, 0.5,
                        0, 0.5;
            ac = pow(2.0*v,1.0/3.0);    
            break;
        case 2:
            unit_cell.resize(3,4);
            unit_cell << 0, 0.5, 0.5,   0,
                        0,   0, 0.5, 0.5,
                        0, 0.5,   0, 0.5;
            ac = pow(4.0*v,1.0/3.0);    
    

            break;
        default:
            throw "Unknown cell type!";
            break;
      }
    }
    catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      exit(1);
    }
    catch(...) {
      cerr << "Exception of unknown type2!\n";
      exit(1);
    } 
    
    CONFIG config(l, ac, grains, v, unit_cell);
    config.name = name;
    config.out_type = out_type;
    config.cell_type = cell_type;

    return config;
}
