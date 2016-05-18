#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <cstdlib>
#include <iomanip>

using namespace std;

struct RES_SIDE
	{
	 int atomn;
	 int atomn_s;
	 int axisn;
	 int *axis;
	 string res_t;
	 string *atomt;
	 string *atomt_s;
	 string **rotate;
	 int *rttn;
	 int **rtt;
	};

class Read_res
	{
	 public:
	 int resn;
	 RES_SIDE *res;
	 int get_res(ifstream *);
	 int get_res(ifstream *,int);
	 int getres(ifstream *);
	 int get_nh();
	 int getside();
	};
