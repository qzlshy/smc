#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

struct Topic
	{
	 string atomid[4];
	 double ic[5];
	};

class Topprm
	{
	 public:
	 int atomn,bondn,imprn,cmapn,icn;
	 string resid;
	 string *atom2;
	 string *atom;
	 double *charge;
	 string **bond;
	 string **impr;
	 string **CMAP;
	 Topic *resic;
	 int getan(ifstream *);
	 int getin(ifstream *);
	 int getbn(ifstream *);
	 int getimn(ifstream *);
	 int getcn(ifstream *);
	 int readatom(ifstream *);
	 int readic(ifstream *);
	 int readbond(ifstream *);
	 int readimpr(ifstream *);
	 int readcmap(ifstream *);
	 int readtop(ifstream *,string);
	 int readtopter(ifstream *,string);
	};
