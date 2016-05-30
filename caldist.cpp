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
 maxr=0.0;
 for(i=0;i<vdw_n;i++)
 if(maxr<vdw_r[i])
	maxr=vdw_r[i];
 goodn=new int[b_s.resnum];
 good=new int*[b_s.resnum];
 for(i=0;i<b_s.resnum;i++)
	good[i]=new int[b_s.rtmn[i]];
 
}

int Caldist::cal(COOR_m &c_m,Cpy_s &cy_s,Rtm_tmp &r_t,int n)
{
 int i,k,t;
 k=0;
 for(i=0;i<r_t.rtmn[n];i++)
	{
	 t=caleach(c_m,cy_s,r_t.cs[n][i],n);
	 if(t==1)
		{
		 good[n][k]=i;
		 k++;
		}
	}
 goodn[n]=k;
}


int Caldist::getw(double a)
{
 w=a;
}

int Caldist::caleach(COOR_m &c_m,Cpy_s &cy_s,COOR_s &cs,int n)
{
 int i,j,k,k1,k2;
 double maxr,cut1,cut2,l1,l2;

 for(i=0;i<c_m.resnum;i++)
 if(i!=n)
	{
	 cut1=maxr*w+maxr*w+c_m.mr[i]+cs.mr;
	 cut1*=cut1;
	 L_L_ab(c_m.cx[i],c_m.cy[i],c_m.cz[i],cs.cx,cs.cy,cs.cz,l1);
	 if(l1<cut1)
		{
		 for(k1=c_m.res_be[i][0];k1<c_m.res_be[i][1];k1++)
			{if(i==n-1)
			 if(c_m.atomt[k1]=="C")
				continue;
			 if(i==n+1)
			 if(c_m.atomt[k1]=="N")
				continue;
			 for(k2=0;k2<cs.atomn;k2++)
				{
				 cut2=vm.r[k1]*w+vs[n].r[k2]*w;
				 cut2*=cut2;
				 L_L_ab(c_m.x[k1],c_m.y[k1],c_m.z[k1],cs.x[k2],cs.y[k2],cs.z[k2],l2);
				 if(l2<cut2)
					return(0);
				}
			}
		}
	}

 for(i=0;i<resnum;i++)
 if(cy_s.ipk[i]==1)
	{
	 cut1=maxr*w+maxr*w+cy_s.cs[i].mr+cs.mr;
	 cut1*=cut1;
	 L_L_ab(cy_s.cs[i].cx,cy_s.cs[i].cy,cy_s.cs[i].cz,cs.cx,cs.cy,cs.cz,l1);
	 if(l1<cut1)
		{
		 for(k1=0;k1<cy_s.cs[i].atomn;k1++)
		 for(k2=0;k2<cs.atomn;k2++)
			{
			 cut2=vs[i].r[k1]*w+vs[n].r[k2]*w;
			 cut2*=cut2;
			 L_L_ab(cy_s.cs[i].x[k1],cy_s.cs[i].y[k1],cy_s.cs[i].z[k1],cs.x[k2],cs.y[k2],cs.z[k2],l2);
			 if(l2<cut2)
				return(0);
			}
		}
	}

 return(1);
 
}
