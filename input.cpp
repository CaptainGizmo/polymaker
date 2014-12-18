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
    string cell, name, out, filename;
    MatrixXd unit_cell;
    Vector3d l;

    try
    {
        // Declare the supported options.
        po::options_description desc("Allowed options");
        desc.add_options()
            ("help", "produce help message")

            ("vol,v",     po::value<double>(&v)->required(),            "volume per atom")
            ("grains,g",  po::value<int>(&grains)->default_value(1),    "number of grains")
            ("width,a",   po::value<double>(&l[0])->required(),         "box side a")
            ("height,b",  po::value<double>(&l[1]),                     "box side b")
            ("length,c",  po::value<double>(&l[2]),                     "box side c")
            ("type,t",    po::value<std::string>(&cell)->default_value("fcc"),                "lattice type")
            ("out,o",     po::value<std::string>(&out)->default_value("lammps"),        "output file type")
            ("input,i",   po::value<std::string>(&filename),   "grain parameters' file")
            ("element,e", po::value<std::string>(&name)->default_value("Xe"),           "element name")
            ("init",      po::value<std::string>()->default_value("Y"),       "")
//            ("time", po::value<int>(), "time")
        ;

        po::variables_map vm;
        po::store(po::parse_command_line(argc, argv, desc), vm);
        po::notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            exit(1);
        }
        if (!vm.count("height")) 
            l[1] = l[0];
        if (!vm.count("length")) 
            l[2] = l[0];
        if (vm.count("out"))     
            cout << "Output file type is " << vm["out"].as<std::string>() << "\n";
          
        if (cell.compare("bcc") == 0)
            cell_type = 1;
        else 
        if (cell.compare("fcc") == 0)
           cell_type = 2;
        else
           throw "Unknown unit cell!";


        if (out.compare("dlpoly") == 0)
            out_type = 1;
        else 
        if (out.compare("lammps") == 0)
           out_type = 2;
        else 
        if (out.compare("xyz") == 0)
           out_type = 3;
        else
           throw "Unknown output format!";
           
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
    
        CONFIG config(l, ac, grains, v, unit_cell);
        config.name = name;
        config.out_type = out_type;
        config.cell_type = cell_type;
        config.cell = cell;
        config.filename = filename;
        if (vm.count("input")) 
            config.read_grains = true;

        if (vm["init"].as<std::string>().compare("Y") == 0)
            config.init = true;
        else 
        if (vm["init"].as<std::string>().compare("N") == 0)
            config.init = false;
        else
           throw "Unknown init param!";


        return config;
    }
    catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      exit(1);
    }
    catch(...) {
      cerr << "Exception of unknown type!\n";
      exit(1);
    } 
}

unsigned ReadGrains(CONFIG &config)
{
    string buf;
    istringstream iss;
    stringstream tmp;
    int i = 0;
    int j;

    try
    {    

        cout << "Read grain parameters from " << config.filename << endl;
        ifstream infile (config.filename.c_str(),ios::in);
        if( !infile ) 
        {
            cout << "Can't open " << config.filename.c_str() << endl; 
            exit(1);
        }

        while( getline( infile, buf) )
        {
            iss.str(buf); 
            iss.clear(); 
            iss >> j >> config.grain[i].r(0)  >> config.grain[i].r(1) >> config.grain[i].r(2) >> config.grain[i].angle(0) >> config.grain[i].angle(1) >> config.grain[i].angle(2);
            //cout << i << " " << config.grain[i].r.transpose() << " " << config.grain[i].angle.transpose() << endl;
            i ++;
        }
        
        if (i != config.grains) throw "Wrong number of lines in grains parameter file.";
        infile.close();

    }
    catch(exception& e) {
      cerr << "error: " << e.what() << "\n";
      exit(1);
    }
/*
    catch(...) {
      cerr << "Exception of unknown type! \n";
      exit(1);
    }
*/
    return 0;
}