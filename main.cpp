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

#ifndef  INPUT_H  
#include "input.h"
#endif

#ifndef  OUTPUT_H  
#include "output.h"
#endif

#ifndef  CLASSES_H  
#include "classes.h"
#endif

using namespace Eigen;
using namespace std;

unsigned IniGrainCenters(CONFIG& config) 
{
    Vector3d a;
    Vector3d rotv;
    Vector3d rtmp;
    double ANGLE;

    cout << "Initialize grain centers" << endl;
    for (int ig=0; ig<config.grains; ig++)
    {
        for (int dim=0; dim<3; dim++)
        {
            config.grain[ig].r(dim) = config.shift(dim) * (rand()%2000000/1000000. - 1);
        }
        for (unsigned j=0; j<3; j++)
    	{
		    config.grain[ig].ang(j) = 2.0 * M_PI * (rand()%1000000/1000000.);
	    }
	    
	    a = config.grain[ig].ang;
	    rotv << 1/cos(a(0))*sin(a(1)), 1/sin(a(0))*sin(a(1)), 1/cos(a(1));
    	ANGLE = -a(2);
	    rtmp = -config.grain[ig].r;
    	config.grain[ig].rotvT = rtmp * cos(ANGLE) + rotv.cross(rtmp) * sin(ANGLE) + (1 - cos(ANGLE)) * rotv *(rotv.dot(rtmp)) + config.grain[ig].r; //  Rodrigue's rotation formula
    }
    return 0;
}


float Lattice(CONFIG &config, int ig)
{
	unsigned ia = 0;
	Vector3d i;
	float N = pow(config.atoms_grain/config.cell.cols(),1.0/3.0);

	for(i(0) = 0; i(0) < N; i(0)++)
	{
		for(i(1) = 0; i(1) < N; i(1)++)
		{
			for(i(2) = 0; i(2) < N; i(2)++)
			{
				for(unsigned col = 0; col < config.cell.cols(); col++)
				    {
					    config.atom_grain[ia+col].r = config.ac * (config.cell.col(col)+i) - config.shift;
					    config.atom_grain[ia+col].type = ig;
					}
				ia+=config.cell.cols();
				//cout << ia << " " << config.atoms_grain << " " << N << " " << i.transpose() <<  endl;
			}
		}
	}
	
	return 0;
}

unsigned RotateBox(CONFIG& config, int ig)
{
	Vector3d a;
	Vector3d rotv;
	Vector3d rtmp;
	double ANGLE;

    a = config.grain[ig].ang;
	
	rotv << cos(a(0))*sin(a(1)), sin(a(0))*sin(a(1)), cos(a(1));
	ANGLE = a(2);
	
	for (int i = 0; i < config.atoms_grain; i++)
	{
	    //cout << i << " " << config.atoms_grain << " " << config.atom_grain.size() <<  endl;
		rtmp = config.atom_grain[i].r;
		config.atom_grain[i].r = rtmp * cos(ANGLE) + rotv.cross(rtmp) * sin(ANGLE) + (1 - cos(ANGLE)) * rotv *(rotv.dot(rtmp)) + config.grain[ig].r; //  Rodrigue's rotation formula
	}

	return 0;
}

unsigned BoxBorder(CONFIG& config, int ig) 
{

        cout << "Box borders check \n";


        Vector3d l0;
        l0 << 0.00001,0.00001,0.00001;
        // check for box borders
        // also here we initialise GoodIndx array

	for(int i = 0; i < config.atoms_grain; i++) 
	{ // check all atoms in grain
		for(int dim=0; dim < 3; dim++)
		{
			if (config.atom_grain[i].r(dim) > config.shift(dim))
				config.atom_grain[i].r(dim) -= 2*config.shift(dim);
			if (config.atom_grain[i].r(dim) < config.shift(dim))
				config.atom_grain[i].r(dim) += 2*config.shift(dim);
		}

	}
	return 0;
}

unsigned Voronoi(CONFIG &config, int ig) 
{

      float r = 0;
      float r0 = 0;
      int jg = 0;
      Vector3d g; g << 0 ,0 , 0;
      Vector3d p; p << 0 ,0 , 0;
      int n = 0;
      float dx = 0.001; //   config.fnn/2 * 0.9;

    #pragma omp parallel for shared(config,n) private(p,g,jg,r,r0)
	for(int i = 0; i < config.atoms_grain; i++) // check all atoms of this grain
	{  
		r = (config.atom_grain[i].r - config.grain[ig].r).norm() + dx; // distance to its center

		for(jg = 0; jg < config.grains; jg++) // compare to the another grain centers
		{
			for(p(0) = -1; p(0) <=1; p(0)++)
			{
				for(p(1) = -1; p(1) <=1; p(1)++)
				{
					for(p(2) = -1; p(2) <=1; p(2)++)
					{
						g = p.array()*config.l.array();
						r0 = (config.atom_grain[i].r - (config.grain[jg].r + g)).norm(); // distance from orig to other center
						if (!((ig == jg) && (p.norm()==0))) // exclude selfcheckong
						if (r > r0) break;
					}
					if(p(2)!=2) break;
				}
				if(p(1)!=2) break;
			}
			if(p(0)!=2) break;
		}
        
        #pragma omp critical
		if (jg==config.grains)
		{
			config.atom_box[config.atoms_box] = config.atom_grain[i];
			config.atoms_box++;
			n++;
		}
		
	}
	cout << "[" << ig << "] " << n << endl;
	
	return 0;
}

unsigned Save(CONFIG &config) 
{

    int dim;

    cout << "Add " << config.atoms_box << " atoms to the box" << endl;

	for(int i = 0; i < config.atoms_box; i++) // check all atoms of this grain
	{
		for(dim=0; dim < 3; dim++)
		{
			if (config.atom_box[i].r(dim) >  config.shift(dim)) config.atom_box[i].r(dim) -= config.l(dim);
			else
			if (config.atom_box[i].r(dim) < -config.shift(dim)) config.atom_box[i].r(dim) += config.l(dim);
		}
	}
	return 0;
}

int main(int argc, char **argv) {

	time_t  time1, time2;
	time(&time1);

	CONFIG config = getinput(argc, argv);

	cout << config.l.transpose() << " box is requested with " << config.ac << " lattice constant and " << config.grains << " grains." << endl;
	//cout << N << " atoms in box, " << config.atoms_grain << " per grain" << endl;
	if (config.init) srand(time(NULL));

	IniGrainCenters(config);
	for (int ig = 0; ig < config.grains; ig++)
	{
		Lattice(config,ig);
		RotateBox(config,ig);
		Voronoi(config,ig);
	}
	cout << endl;
	Save(config);
	
	printLAMMPS(config);
	printDLPOLY(config);
	printXYZ(config);
	
	time(&time2);
	
	cout << "Run " << time2-time1 << " s." << endl;

	return 0;
}
