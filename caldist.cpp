#include "long.h"

string vdw_a[5]={"H","C","N","O","S"};
double vdw_r[5]={1.2,1.7,1.55,1.52,1.85};
int vdw_n=5;

int Caldist::initm(Build_m &b_m)
{
 int i,j;
 string tmp_v;
 vm.atomn=b_m.atomn;
 vm.r=new double[vm.atomn];
 for(i=0;i<vm.atomn;i++)
	{tmp_v=b_m.atomt[i].substr(0,1);
	 for(j=0;j<vdw_n;j++)
	 if(tmp_v==vdw_a[j])
		{
		 vm.r[i]=vdw_r[j];
		 break;
		}
	}
}

int Caldist::inits(Build_s &b_s)
{
 int i,j,k;
 string tmp_v;
 resnum=b_s.resnum;
 vs=new VDW_r[resnum];
 for(i=0;i<resnum;i++)
	{
	 vs[i].atomn=b_s.ss[i].atomn;
	 vs[i].r=new double[vs[i].atomn];
	 for(j=0;j<vs[i].atomn;j++)
		{tmp_v=b_s.ss[i].atomt[j].substr(0,1);
		 for(k=0;k<vdw_n;k++)
		 if(tmp_v==vdw_a[k])
			{
			 vs[i].r[j]=vdw_r[k];
			 break;
			}
		}
	}
}

int Caldist::init(Build_m &b_m,Build_s &b_s)
{
 int i,j;
 initm(b_m);
 inits(b_s);
 goodn=new int[b_s.resnum];
 good=new int*[b_s.resnum];
 for(i=0;i<b_s.resnum;i++)
	good[i]=new int[b_s.rtmn[i]];
 
}

int Caldist::cal(COOR_m &c_m,Cpy_s &cy_s,Rtm_tmp &r_t,int n)
{
 int i,k,t;
 k=0;
 for(i=0;i<rtmn[n];i++)
	{
	 t=caleach(c_m,cy_s,r_t.cs[n][i]);
	 if(t==1)
		{
		 good[n][k]=i;
		 k++;
		}
	}
 goodn[n]=k;
}

int Caldist::caleach(COOR_m &c_m,Cpy_s &cy_s,COOR_s &cs)
{
 
}
