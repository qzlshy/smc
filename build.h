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
	 int atomn;
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
	 int atomn;
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
	 int build(Build_m,Read_res,Res_rtm);
	 int get_ss(Build_m,Read_res);
	 int get_cs(Res_rtm);
	 int rt_fit(Build_m);
	 int fit_one(Build_m,int,int);
	 int resnum;
	 int *rtmn;
	 string *res_t;
	 Structure_s *ss;
	 COOR_s **cs;
	};
