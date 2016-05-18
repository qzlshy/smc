#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>
#include <cstdlib>
#define _USE_MATH_DEFINES
#include <cmath>

using namespace std;

#define L_CHAX(ax,ay,az,bx,by,bz,cx,cy,cz) { \
        cx=(ay)*(bz)-(by)*(az); \
        cy=(az)*(bx)-(bz)*(ax); \
        cz=(ax)*(by)-(bx)*(ay);            }

#define L_DIANX(ax,ay,az,bx,by,bz,c)      { \
        c=(ax)*(bx)+(ay)*(by)+(az)*(bz);  }

#define L_MO(ax,ay,az,b)                       { \
        b=sqrt((ax)*(ax)+(ay)*(ay)+(az)*(az)); }

#define L_MO_ab(ax,ay,az,bx,by,bz,b)                       { \
        b=sqrt((ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz)); }

#define L_L_ab(ax,ay,az,bx,by,bz,b)                       { \
        b=(ax-bx)*(ax-bx)+(ay-by)*(ay-by)+(az-bz)*(az-bz); }

int rotation_a(double *,double *,double *,int,double *,double);
int rotation_n(double *,double *,double *,int,int *,double *,double);
int rotation(double *,double *,double);
double get_angle(double,double,double,double,double,double);
double get_dhd(double *,double *,double *);
class Guess_d
        {
         public:
         int init(double *,double *,double *,double,double,double);
         int rotation(double *,double *,double);
         int guess();
         int getd(double *);
         private:
         double a[3];
         double b[3];
         double c[3];
         double d[3];
         double cd,bcd,abcd;
        };
