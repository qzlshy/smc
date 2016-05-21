#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

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

struct COOR_m
	{
	 public:
	 int init_cm(Build_m &);
	 int atomn,resnum;
	 double *x,*y,*z;
	 double **N,**CA,**C,**CB;
	 double *cx,*cy,*cz,*mr;
	};

class COOR_s
	{
	 public:
	 int atomn;
	 int operator=(COOR_s &);
	 double *x,*y,*z;
	 double *N,*CA,*C,*CB;
	 double cx,cy,cz,mr;
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
	 int get_ct(int,int);
	 int resnum;
	 int *rtmn;
	 string *res_t;
	 Structure_s *ss;
	 COOR_s **cs;
	};

class Cpy_s
	{
	 public:
	 int resnum;
	 int ipk;
	 COOR_s *cs;
	 int init_cs(Build_s &);
	};
