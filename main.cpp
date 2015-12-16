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
    double angle;
    time_t  time1, time2;
    time(&time1);


    cout << "Initialize grain centers ";
    for (int ig=0; ig<config.grains; ig++)
    {
        for (int dim=0; dim<3; dim++)
        {
            config.grain[ig].r(dim) = config.shift(dim) * (rand()%2000000/1000000. - 1);
        }

        for (unsigned j=0; j<3; j++)
        {
            config.grain[ig].angle(j) = 2.0 * M_PI * (rand()%1000000/1000000.);
        }

        a = config.grain[ig].angle;
        rotv << 1/cos(a(0))*sin(a(1)), 1/sin(a(0))*sin(a(1)), 1/cos(a(1));
        angle = -a(2);
        rtmp = -config.grain[ig].r;
        config.grain[ig].rotvT = rtmp * cos(angle) + rotv.cross(rtmp) * sin(angle) + (1 - cos(angle)) * rotv *(rotv.dot(rtmp)) + config.grain[ig].r; //  Rodrigue's rotation formula
        cout << ".";
    }

    time(&time2);
    if (config.time) cout << endl << "Done in " << time2-time1 << " s.";
    cout << endl;
    return 0;
}


float Lattice(CONFIG &config, int ig)
{
    unsigned ia = 0;
    Vector3d i;
    time_t  time1, time2;
    time(&time1);

    float N = pow(config.atoms_grain/config.unit_cell.cols(),1.0/3.0);

    for(i(0) = 0; i(0) < N; i(0)++)
    {
        for(i(1) = 0; i(1) < N; i(1)++)
        {
            for(i(2) = 0; i(2) < N; i(2)++)
            {
                for(unsigned col = 0; col < config.unit_cell.cols(); col++)
                {
                    config.atom_grain[ia+col].r = config.ac * (config.unit_cell.col(col)+i) - config.shift;
                    config.atom_grain[ia+col].type = ig;
                }
                ia+=config.unit_cell.cols();
            }
        }
    }

    time(&time2);
    if (config.time) cout << endl << "Lattice init " << time2-time1 << " s. ";
    return 0;
}

unsigned RotateBox(CONFIG& config, int ig)
{
    Vector3d a;
    Vector3d rotv;
    Vector3d rtmp;
    double   angle;
    time_t   time1, time2;
    time(&time1);


    a = config.grain[ig].angle;

    rotv << cos(a(0))*sin(a(1)), sin(a(0))*sin(a(1)), cos(a(1));
    angle = a(2);

    for (int i = 0; i < config.atoms_grain; i++)
    {
        rtmp = config.atom_grain[i].r;
        config.atom_grain[i].r = rtmp * cos(angle) + rotv.cross(rtmp) * sin(angle) + (1 - cos(angle)) * rotv *(rotv.dot(rtmp)) + config.grain[ig].r; //  Rodrigue's rotation formula
    }

    time(&time2);
    if (config.time) cout << "Rotated in " << time2-time1 << " s.";

    return 0;
}

unsigned Voronoi(CONFIG &config, int ig)
{

    float r  = 0;
    float r0 = 0;
    int   jg = 0;
    Vector3d g; g << 0 ,0 , 0;
    Vector3d p; p << 0 ,0 , 0;
    int n = 0;
    time_t  time1, time2;
    time(&time1);

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

    cout << "[" << ig << "] " << n ;
    time(&time2);
    cout << " Done in " << time2-time1 << " s.";
    cout << endl;

    return 0;
}

unsigned Save(CONFIG &config)
{
    time_t  time1, time2;
    time(&time1);
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

    time(&time2);
    if (config.time) cout << "Done in " << time2-time1 << " s." << endl;

    return 0;
}

unsigned SortGrig(CONFIG &config)
{
    Vector3d Lcell;
    Lcell = config.l / config.grid_size;

    for(int i=0; i < config.atoms_box; i++)
    {
        // shift coords to positive values and divide by grid cell length
        // Eigen::ceil() appears from v3.3
        Vector3d r_reduced = (config.atom_box[i].r + config.shift).array() / Lcell.array();
        // convert reduced coords to linear array number
        int grid_id = floor(r_reduced[0]) * config.grid_size * config.grid_size + floor(r_reduced[1]) * config.grid_size + floor(r_reduced[2]);
        //cout << grid_id <<" "<< floor(r_reduced[0]) << " " << floor(r_reduced[1]) << " " << floor(r_reduced[2]) << endl;
        config.grid[grid_id].id.push_back(i);
    }

    return 0;
}

int main(int argc, char **argv) {

    time_t  time1, time2;
    time(&time1);

    CONFIG config = getinput(argc, argv);

    cout << config.l.transpose() << " box, " << config.cell << " is requested, " << config.v << " A^3 per atom (lattice constant " << config.ac << ") in " << config.grains << " grains." << endl;
    //cout << N << " atoms in box, " << config.atoms_grain << " per grain" << endl;

    if (config.read_grains)
    {
        ReadGrains(config);
        srand(time(NULL));
    }
    else
    {
        if (config.init) srand(time(NULL));
        IniGrainCenters(config);
        WriteGrains(config);
    }
    for (int ig = 0; ig < config.grains; ig++)
    {
        Lattice(config,ig);
        RotateBox(config,ig);
        Voronoi(config,ig);
    }
    cout << endl;
    Save(config);
    SortGrig(config);

    switch(config.out_type)
    {
        case 1:
            WriteDLPOLY(config);
            break;
        case 2:
            WriteLAMMPS(config);
            break;
        case 3:
            WriteXYZ(config);
            break;
        default:
            throw "Unknown output format.";
            break;
    }

    time(&time2);

    cout << "Total run " << time2-time1 << " s." << endl;

    return 0;
}
