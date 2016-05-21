#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

class VDW_r
	{
	 public:
	 int atomn;
	 double *r;
	};

class Caldist
	{
	 public:
	 int initm(Build_m &);
	 int inits(Build_s &);
	 int init(Build_m &,Build_s &);
	 int resnum;
	 VDW_r vm;
	 VDW_r *vs;
	};
