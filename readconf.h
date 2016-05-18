#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

class Read_conf
	{
	 public:
	 Read_conf(char *);
	 int inputconf(char *);
	 string getstr(char *);
	 int getint(char *);
	 int getfloat(char *);
	 private:
	 int n;
	 string *line;
	};
