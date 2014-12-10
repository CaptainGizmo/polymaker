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

    int grains, N;
    double ac, v;
    double a, b, c;

    try
    {
      // Declare the supported options.
      po::options_description desc("Allowed options");
      desc.add_options()
          ("help", "produce help message")

          //("optimization",   po::value<int>(&opt)->default_value(10),    "optimization level")
          //("include-path,I", po::value< vector<string> >()->composing(), "include path")

          ("V,v", po::value<double>(), "volume per atom")
          ("N,n", po::value<int>(), "number of grains")
          ("A,a", po::value<double>()->composing(), "box side a")
          ("B,b", po::value<double>(), "box side b")
          ("C,c", po::value<double>(), "box side c")
          ("type,t", po::value< vector<string> >(), "lattice type")
          ("out,o", po::value< vector<string> >(), "output file type")
      ;
      po::variables_map vm;
      po::store(po::parse_command_line(argc, argv, desc), vm);
      po::notify(vm);    

      if (vm.count("help")) {
        cout << desc << "\n";
        exit(1);
      }
      
      if (vm.count("out"))
      {
          cout << "Output file type is " << vm["out"].as< vector<string> >() << "\n";
      }
    
      if (vm.count("A")) 
      {
          a = vm["A"].as<double>();
      } 
      else
      {
          cout << "Please set box size a.\n";
          exit(1);
      }     

      if (vm.count("B")) 
      {
          b = vm["B"].as<double>();
      } 
      else
      {
          b = a;
      }

      if (vm.count("C")) 
      {
          c = vm["C"].as<double>();
      } 
      else
      {
          c = a;
      }
    
      if (vm.count("N")) 
      {
          grains = vm["N"].as<int>();
      }    
      else
      {
          grains = 1;
      }


      if (vm.count("V")) 
      {
          v = vm["V"].as<double>();
      } 
      else
      {
          v = 10;
      }
    
    }
    catch(exception& e) {
        cerr << "error: " << e.what() << "\n";
        exit(1);
    }
    catch(...) {
        cerr << "Exception of unknown type!\n";
        exit(1);
    } 
    
    Vector3d l(a,b,c);
    ac = pow(4.0*v,1.0/3.0);      // bcc
    N = ceil(l.prod()/v);
    CONFIG config(l,ac,grains,N);
    
    return config;
}