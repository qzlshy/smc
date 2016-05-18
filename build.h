#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

struct COOR_m
	{
	 double *x,*y,*z;
	 double **N,**CA,**C,**CB;
	};

class Build_m
	{
	 public:
	 int build(Main_chain,Topprm *);
	 int get_mc(Main_chain);
	 int guess_CB_miss(Topprm *);
	 int atomn,resnum;
	 int *res_id;
	 string *atomt,*res_t;
	 char *res_chain;
	 double *x,*y,*z;
	 int **res_be;
	 int *CB_miss;
	 double **N,**CA,**C,**CB;
	};

class COOR_s
	{
	 public:
	 double *x,*y,*z;
	 double N[3],CA[3],C[3],CB[3];
	};

class Structure_s
	{
	 public:
	 int atomn;
	 string *atomt;
	 int axisn;
	 int *axis;
	 int *rttn;
	 int **rtt;
	};

class Build_s
	{
	 public:
	 int get_ss(Main_chain,Read_res);
	 int resnum;
	 string *res_t;
	 Structure_s *ss;
	};
