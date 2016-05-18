#include "long.h"

int rotation_a(double *r0,double *r1,double *r2,int n,double *v,double theta)
{
 int i;
 double t1[3],t2[3],t3[3];
 double vr;
 double c_theta,s_theta;
 c_theta=cos(theta);
 s_theta=sin(theta);

 for(i=0;i<n;i++)
        {
         t1[0]=r0[i]*c_theta; t1[1]=r1[i]*c_theta; t1[2]=r2[i]*c_theta;
         L_CHAX(v[0],v[1],v[2],r0[i],r1[i],r2[i],t2[0],t2[1],t2[2]);
         t2[0]=t2[0]*s_theta; t2[1]=t2[1]*s_theta; t2[2]=t2[2]*s_theta;
         L_DIANX(v[0],v[1],v[2],r0[i],r1[i],r2[i],vr);
         t3[0]=vr*v[0]*(1-c_theta); t3[1]=vr*v[1]*(1-c_theta); t3[2]=vr*v[2]*(1-c_theta);

         r0[i]=t1[0]+t2[0]+t3[0];
         r1[i]=t1[1]+t2[1]+t3[1];
         r2[i]=t1[2]+t2[2]+t3[2];
        }

}

int rotation_n(double *r0,double *r1,double *r2,int n1,int *n2,double *v,double theta)
{
 int i;
 double t1[3],t2[3],t3[3];
 double vr;
 double c_theta,s_theta;
 c_theta=cos(theta);
 s_theta=sin(theta);

 for(i=0;i<n1;i++)
        {
         t1[0]=r0[n2[i]]*c_theta; t1[1]=r1[n2[i]]*c_theta; t1[2]=r2[n2[i]]*c_theta;
         L_CHAX(v[0],v[1],v[2],r0[n2[i]],r1[n2[i]],r2[n2[i]],t2[0],t2[1],t2[2]);
         t2[0]=t2[0]*s_theta; t2[1]=t2[1]*s_theta; t2[2]=t2[2]*s_theta;
         L_DIANX(v[0],v[1],v[2],r0[n2[i]],r1[n2[i]],r2[n2[i]],vr);
         t3[0]=vr*v[0]*(1-c_theta); t3[1]=vr*v[1]*(1-c_theta); t3[2]=vr*v[2]*(1-c_theta);

         r0[n2[i]]=t1[0]+t2[0]+t3[0];
         r1[n2[i]]=t1[1]+t2[1]+t3[1];
         r2[n2[i]]=t1[2]+t2[2]+t3[2];
        }

}

int rotation(double *r,double *v,double theta)
{
 double t1[3],t2[3],t3[3];
 double vr;
 t1[0]=r[0]*cos(theta); t1[1]=r[1]*cos(theta); t1[2]=r[2]*cos(theta);
 L_CHAX(v[0],v[1],v[2],r[0],r[1],r[2],t2[0],t2[1],t2[2]);
 t2[0]=t2[0]*sin(theta); t2[1]=t2[1]*sin(theta);t2[2]=t2[2]*sin(theta);
 L_DIANX(v[0],v[1],v[2],r[0],r[1],r[2],vr);
 t3[0]=vr*v[0]*(1-cos(theta)); t3[1]=vr*v[1]*(1-cos(theta)); t3[2]=vr*v[2]*(1-cos(theta));

 r[0]=t1[0]+t2[0]+t3[0];
 r[1]=t1[1]+t2[1]+t3[1];
 r[2]=t1[2]+t2[2]+t3[2];
}

double get_angle(double ax,double ay,double az,double bx,double by,double bz)
{
 double a,b,c,tmp;
 L_DIANX(ax,ay,az,bx,by,bz,c);
 L_MO(ax,ay,az,a);
 L_MO(bx,by,bz,b);
 tmp=c/(a*b);
 if(tmp>1.0) tmp=1.0;
 if(tmp<-1.0) tmp=-1.0;
 return(acos(tmp));
}

double get_dhd(double *a,double *b,double *c)
{
 double fab[3],fbc[3],ff[3];
 double m_fab,m_fbc;
 double angle,tmp,d;

 L_CHAX(a[0],a[1],a[2],b[0],b[1],b[2],fab[0],fab[1],fab[2]);
 L_CHAX(b[0],b[1],b[2],c[0],c[1],c[2],fbc[0],fbc[1],fbc[2]);
 L_MO(fab[0],fab[1],fab[2],m_fab);
 L_MO(fbc[0],fbc[1],fbc[2],m_fbc);
 L_DIANX(fab[0],fab[1],fab[2],fbc[0],fbc[1],fbc[2],d);
 tmp=d/(m_fab*m_fbc);
 if(tmp>1.0) tmp=1.0;
 if(tmp<-1.0) tmp=-1.0;
 angle=acos(tmp);

 L_CHAX(fab[0],fab[1],fab[2],fbc[0],fbc[1],fbc[2],ff[0],ff[1],ff[2]);
 if((ff[0]*b[0]+ff[1]*b[1]+ff[2]*b[2])<0)
 angle=-angle;

 return(angle);
}

int Guess_d::init(double *ma,double *mb,double *mc,double mabcd,double mbcd,double mcd)
{
 a[0]=ma[0]; a[1]=ma[1]; a[2]=ma[2];
 b[0]=mb[0]; b[1]=mb[1]; b[2]=mb[2];
 c[0]=mc[0]; c[1]=mc[1]; c[2]=mc[2];
 cd=mcd; bcd=mbcd*M_PI/180.0; abcd=mabcd*M_PI/180.0;
}

int Guess_d::rotation(double *r,double *v,double theta)
{
 double t1[3],t2[3],t3[3];
 double vr;
 t1[0]=r[0]*cos(theta); t1[1]=r[1]*cos(theta); t1[2]=r[2]*cos(theta);
 L_CHAX(v[0],v[1],v[2],r[0],r[1],r[2],t2[0],t2[1],t2[2]);
 t2[0]=t2[0]*sin(theta); t2[1]=t2[1]*sin(theta);t2[2]=t2[2]*sin(theta);
 L_DIANX(v[0],v[1],v[2],r[0],r[1],r[2],vr);
 t3[0]=vr*v[0]*(1-cos(theta)); t3[1]=vr*v[1]*(1-cos(theta)); t3[2]=vr*v[2]*(1-cos(theta));

 r[0]=t1[0]+t2[0]+t3[0];
 r[1]=t1[1]+t2[1]+t3[1];
 r[2]=t1[2]+t2[2]+t3[2];
}

int Guess_d::guess()
{
 double ab[3],cb[3],bc[3],fabc[3];
 double fmo,cbmo,bcmo;

 ab[0]=b[0]-a[0];
 ab[1]=b[1]-a[1];
 ab[2]=b[2]-a[2];
 cb[0]=b[0]-c[0];
 cb[1]=b[1]-c[1];
 cb[2]=b[2]-c[2];
 bc[0]=c[0]-b[0];
 bc[1]=c[1]-b[1];
 bc[2]=c[2]-b[2];

 L_MO(cb[0],cb[1],cb[2],cbmo);
 d[0]=cb[0]*cd/cbmo; d[1]=cb[1]*cd/cbmo; d[2]=cb[2]*cd/cbmo;

 L_MO(bc[0],cb[1],bc[2],bcmo);
 bc[0]/=bcmo; bc[1]/=bcmo; bc[2]/=bcmo;

 L_CHAX(ab[0],ab[1],ab[2],cb[0],cb[1],cb[2],fabc[0],fabc[1],fabc[2]);
 L_MO(fabc[0],fabc[1],fabc[2],fmo);
 fabc[0]/=fmo; fabc[1]/=fmo; fabc[2]/=fmo;

 rotation(d,fabc,bcd);
 rotation(d,bc,abcd);

 d[0]+=c[0];
 d[1]+=c[1];
 d[2]+=c[2];

}

int Guess_d::getd(double *md)
{
 md[0]=d[0];
 md[1]=d[1];
 md[2]=d[2];
}

