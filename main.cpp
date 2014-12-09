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

unsigned IniGrainCenters(CONFIG& config) 
{
    cout << "Initialize grain centers" << endl;
    for (int ig=0; ig<config.grains; ig++)
    {
        for (int dim=0; dim<3; dim++)
        {
            config.grain[ig].r(dim) = config.shift(dim) * (rand()%2000000/1000000. - 1);
        }
//        cout << config.grain[ig].r.transpose() << endl;
    }
    return 0;
}

float bcc(CONFIG &config)
{

	// generate bcc coordinates, 2 atoms per cell
	Matrix<double, 3, 2> cell;
	
	cell << 0, 0.5,
	        0, 0.5,
	        0, 0.5;

//	cout << "Generating atoms' coordinates" << endl;
        //#pragma omp parallel for
	unsigned ia = 0;
	Vector3d i;
	float N = pow(config.atoms_grain/2.0,1.0/3.0);

	for(i(0) = 0; i(0) < N; i(0)++)
	{
		for(i(1) = 0; i(1) < N; i(1)++)
		{
			for(i(2) = 0; i(2) < N; i(2)++)
			{
				for(unsigned col = 0; col < cell.cols(); col++)
					config.atom_grain[ia+col].r = config.ac * (cell.col(col)+i) - config.shift;
				ia+=2;
			}
		}
	}
	

	return 0;
}



float fcc(CONFIG &config)
{

	// generate fcc coordinates, 4 atoms per cell
	Matrix<double, 3, 4> cell;
	
	cell << 0, 0.5, 0.5,   0,
	        0,   0, 0.5, 0.5,
	        0, 0.5,   0, 0.5; 

//	cout << "Generating atoms' coordinates" << endl;
        //#pragma omp parallel for
	unsigned ia = 0;
	Vector3d i;
	float N = pow(config.atoms_grain/4.0,1.0/3.0);

	for(i(0) = 0; i(0) < N; i(0)++)
	{
		for(i(1) = 0; i(1) < N; i(1)++)
		{
			for(i(2) = 0; i(2) < N; i(2)++)
			{
				for(unsigned col = 0; col < cell.cols(); col++)
					config.atom_grain[ia+col].r = config.ac * (cell.col(col)+i) - config.shift;
				ia+=4;
			}
		}
	}
	

	return 0;
}


unsigned RotateBox(CONFIG& config, int ig)
{
	Vector3d rotv;
	Vector3d rtmp;

//    cout << "Rotate grains" << endl;
	for (unsigned j=0; j<3; j++)
	{
		rotv(j) = -1 + 2*rand()%1000000/1000000.;
	}
	rotv /= rotv.norm();
	
	//rotv << 0, 1, 0;

	 // Generate random angle
	float ANGLE = 2.0 * M_PI * (rand()%1000000/1000000.);

	for (int i = 0; i<config.atoms_grain ;i++)
	{
		rtmp = config.atom_grain[i].r;
		config.atom_grain[i].r = rtmp * cos(ANGLE) + rotv.cross(rtmp) * sin(ANGLE) + (1 - cos(ANGLE)) * rotv *(rotv.dot(rtmp)) + config.grain[ig].r; //  Rodrigues' rotation formula
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

      float r,r0;
      int jg;
      Vector3d g,p;
      int n=0;
      float dx = 1.5; //config.fnn / 2 * 0.9;

//        cout << "Voronoi check \n";

	for(int i = 0; i < config.atoms_grain; i++) // check all atoms of this grain
	{
		r = (config.atom_grain[i].r - config.grain[ig].r).norm() + dx; // distance to its center
		for(jg = 0; jg < config.grains; jg++) // compare to the another grain centers
		{
			//if (ig != jg) // do not need to check to itself
			//{
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
			//}
		}

		if (jg==config.grains)
		{
			config.atom_box[config.atoms_box].r = config.atom_grain[i].r;
			config.atoms_box++;
			n++;
		}
	}
	cout << "[" << ig << "] " << n << endl;
	
	return 0;
}
/*
unsigned DelVoid(CONFIG &config)
{


	Matrix<double, 3, 4> cell;
	
	cell << 0, 0.5, 0.5,   0,
	        0,   0, 0.5, 0.5,
	        0, 0.5,   0, 0.5; 
	        

	unsigned ia = 0;
	Vector3d i;
	Vector3d N;
	N = config.l/config.ac; // not intejer but this is OK!
	vector<ATOM> fill(ceil(4*N.prod()));

	for(i(0) = 0; i(0) < N(0); i(0)++)
	{
		for(i(1) = 0; i(1) < N(1); i(1)++)
		{
			for(i(2) = 0; i(2) < N(2); i(2)++)
			{
				for(unsigned col = 0; col < cell.cols(); col++)
					fill[ia+col].r = config.ac * (cell.col(col)+i)-config.shift;
				ia+=4;
			}
		}
	}

	double r,d;
	d = 0.9 * (config.ac * pow(3,0.5)/2.0);
	int N0 = config.atoms_box;
	int j,k;

	cout << "Filling empty spaces, ";

	for (k = 0; k < 4*N.prod() ; k++)
	{
		for (j = 0; j < N0; j++)
		{
			r = (config.atom_box[j].r - fill[k].r).norm();
			if (r<=d) break; 
		}
		if (j==N0)
		{
			config.atom_box[config.atoms_box].r = fill[k].r;
			config.atoms_box++;
		}
	}
	
	cout << config.atoms_box - N0 << " added." << endl;

	return 0;
}

unsigned Overlap(CONFIG &config)
{

	vector<ATOM> tmp=config.atom_box;
	int i, j, N0,k;
	double r,d;
	d = 0.9 * (config.ac * pow(2,0.5)/2.0); //fcc
	//d = 0.9 * (config.ac * pow(3,0.5)/2.0); //bcc
	N0 = config.atoms_box;
	k=0;

	cout << "Check of overlapped atoms, distance " << d << " A, ";

	for (i = 0 ; i < config.atoms_box; i++)
	{
		for (j = i+1; j < config.atoms_box; j++)
		{
			r = (config.atom_box[j].r - config.atom_box[i].r).norm();
			if (r<d) break;
		}
		if (j==config.atoms_box)
		{
			tmp[k]=config.atom_box[i];
			k++;
		}
//		else
//		{
//			cout << i << " " << j << endl;
//		}
	}

	cout << config.atoms_box - k << " deleted." << endl;
	config.atoms_box=k;
	config.atom_box=tmp;
	return 0;
}

*/

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

unsigned printLAMMPS(CONFIG &config)
{
    cout << "Saving LAMMPS config with " << config.atoms_box << " atoms." << endl;
    stringstream filename;
    string filename_s;
    

    filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains << ".in";
    filename_s = filename.str();
    ofstream outfile (filename_s.c_str(), ios::trunc);
    outfile.precision(5);
    
    Vector3d l;
    
    //outfile << "Lattice file generated by lmpscreator utility for "<< pow(config.ac,3)/4.0 <<" per atom\n\n";
    outfile << "Lattice file generated by lmpscreator utility for "<< pow(config.ac,3)/2.0 <<" per atom\n\n";
    outfile << "1 atom types\n";
    outfile << config.atoms_box << " atoms\n\n";
    outfile << -config.shift(0) << " " << config.shift(0) << " xlo xhi\n";
    outfile << -config.shift(1) << " " << config.shift(1) << " ylo yhi\n";
    outfile << -config.shift(2) << " " << config.shift(2) << " zlo zhi\n\n";
    //outfile << "Pair Coeffs\n\n";
    //outfile << "2 0.42276 4.099\n";
    outfile << endl;
    outfile << "Masses\n\n";
    outfile << "1 55.845\n";
    //outfile << "2 60\n";
    //outfile << "3 131.2937\n";
    outfile << endl;
    outfile << "Atoms\n\n";

    //pair_style      buck 8.0
    //pair_coeff      1 1 177079.254302103 0.343846  6917.93499044

    for (int i = 0; i < config.atoms_box ; i++)
    {
        outfile << i+1 << " 1 " << setw (15)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << endl;
    }

    outfile.close();

    return 0;
}


unsigned printXYZ(CONFIG &config)
{
	stringstream filename;
	string filename_s;

	filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains   << ".xyz";
	filename_s = filename.str();
	ofstream outfile (filename_s.c_str(), ios::trunc);
	outfile.precision(5);

	outfile << config.atoms_box << "\n Fenon\n";
	
	
	for (int  i = 0; i < config.atoms_box; i++) 
	{
		outfile << "Fe" << setw (15)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << endl;
	}

	outfile.close();
	
	return 0;
}

unsigned printDLPOLY(CONFIG &config)
{
	stringstream filename;
	string filename_s;

	filename << "CONFIG." << config.l(0) << "_" << config.l(1) << "_" << config.l(2) << "_" << config.grains   << "";
	filename_s = filename.str();
	ofstream outfile (filename_s.c_str(), ios::trunc);
	outfile.precision(5);

	outfile << "Lattice file generated by lmpscreator utility for "<< pow(config.ac,3)/2.0 <<" per atom" << endl;
	outfile << "0         3" << endl;
	outfile << config.l(0) << " 0.0000000000 " << " 0.0000000000"  <<endl;
	outfile << "0.0000000000 " << config.l(1) << " 0.0000000000"  <<endl;
	outfile << "0.0000000000 " << " 0.0000000000 "  << config.l(2) <<endl;
	
	
	for (int  i = 0; i < config.atoms_box; i++) 
	{
		outfile << "Fe   " << i+1 << endl;
		outfile << setw (15)  << config.atom_box[i].r(0) << " " << config.atom_box[i].r(1) << " " << config.atom_box[i].r(2) << endl;
	}

	outfile.close();
	
	return 0;
}


int main(int argc, char *argv[]) {

	time_t  time1, time2;
	time(&time1);


	Vector3d l(100,100,100);
	int grains=3;
	double v = 10; // cubic A per at

	if (argc==1) {cout << "a b c v n" << endl; exit(0);}
	if (argc>1) l(0) = atof(argv[1]);
	if (argc>2) l(1) = atof(argv[2]);
	if (argc>3) l(2) = atof(argv[3]);
	if (argc>4) v = atof(argv[4]);
	if (argc>5) grains = atoi(argv[5]);

	double ac = pow(4.0*v,1.0/3.0);      // bcc


	int N = ceil(l.prod()/v);

	if ( grains < 10 ) {N *= 32;}
/*
	else
	if ( grains > 1000 ) {N /= 128; }
	else
	if ( grains > 100 ) {N /= 32; }
	else
	if ( grains > 50 ) {N /= 8;}
	N=ceil(N);
*/
	CONFIG config(l,ac,grains,N);
	//config.fnn = ac * pow(3,0.5)/2.0;    // bcc

	cout << config.l.transpose() << " box is requested with " << config.ac << " lattice constant and " << config.grains << " grains." << endl;
	cout << N << " atoms in box, " << config.atoms_grain << " per grain" << endl;
	//srand(time(NULL));

	IniGrainCenters(config);
	for (int ig = 0; ig < grains; ig++)
	{
		fcc(config);
		RotateBox(config,ig);
		//BoxBorder(config);
		Voronoi(config,ig);
//		cout << ".";
	}
	cout << endl;
	Save(config);
	
//	DelVoid(config);
//	Overlap(config);

	printLAMMPS(config);
	printDLPOLY(config);
	printXYZ(config);
	
	time(&time2);
	
	cout << "Run " << time2-time1 << " s." << endl;

	return 0;
}
