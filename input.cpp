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

#include <boost/program_options.hpp>

#ifndef  CLASSES_H  
#include "classes.h"
#endif

using namespace Eigen;
using namespace std;
namespace po = boost::program_options;

unsigned getinput(int argc, const char **argv, CONFIG &config)
{

    // Declare the supported options.
    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("v", po::value<double>(), "volume per atom")
        ("n", po::value<int>(), "number of grains")
        ("a", po::value<double>(), "box side")
    ;
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);    

    if (vm.count("help")) {
        cout << desc << "\n";
        return 1;
    }

    if (vm.count("n")) 
    {
        cout << "Number of grains was set to "  << vm["n"].as<int>() << ".\n";
    } 
    else 
    {
        cout << "Number of grains is set to 1 (default).\n";
    }
    
    return 0;
}
