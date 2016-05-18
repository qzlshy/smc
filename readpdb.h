#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

class Read_pdb
        {
	 public:
	 int readpdb(ifstream *);
	 int getprotein();
	 int choosealt();
	 int removeh();
	 int atomn;
	 string *atomt,*residt,*segid;
	 char *chain;
	 char *alt;
	 int *serial,*resid,*index;
	 double *x,*y,*z;
	 double *tof;
	 string *line;
	};

class Pdb_res: public Read_pdb
	{
	 public:
	 int getres();
	 int resnum;
	 int **res_be;
	 string *res_t;
	 char *res_chain;
	};

struct COOR
        {
	 int n;
	 string *atomt;
         double *x;
         double *y;
         double *z;
         double N[3],CA[3],C[3],CB[3];
        };


class Res_rtm: private Pdb_res
	{
	 public:
	 int get_rtm(ifstream *);
	 int rtn;
	 int *rtmn;
	 string *rtmrt;
	 COOR **rtmcr;
	};

class Main_chain: public Pdb_res
	{
	 public:
	 int get_mc(ifstream *);
	 int get_mc(ifstream *,int);
	 double **N,**CA,**C,**CB;
	 int *N_n,*CA_n,*C_n,*CB_miss;
	};
