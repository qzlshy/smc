#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

struct BONDS
{
 string a1;
 string a2;
 double Kb;
 double b0;
};

struct ANGLES
{
 string a1;
 string a2;
 string a3;
 double Ktheta;
 double Theta0;
 double Kub;
 double S0;
};

struct DIHEDRALS
{
 string a1;
 string a2;
 string a3;
 string a4;
 double Kchi;
 int n;
 double delta;
};

struct IMPROPER
{
 string a1;
 string a2;
 string a3;
 string a4;
 double Kpsi;
 int n;
 double psi0;
};

struct NONBONDED
{
 string a1;
 double ignored1;
 double epsilon;
 double Rmin2;
 double ignored2;
 double eps14;
 double Rmin214;
};

class Read_prm
	{
	 public:
	 int nbond,nangle,ndihedral,nimpro,nnonbond;
	 BONDS *bond;
	 ANGLES *angle;
	 DIHEDRALS *dihedral;
	 IMPROPER *improp;
	 NONBONDED *nonbond;
	 int get_prm(ifstream *);
	};

