#include "long.h"

int Build_m::get_mc(Main_chain m_c)
{
 int i,j,k;
 atomn=m_c.atomn;
 resnum=m_c.resnum;
 atomt=new string[atomn];
 res_id=new int[resnum];
 res_t=new string[resnum];
 res_chain=new char[resnum];
 x=new double[atomn];
 y=new double[atomn];
 z=new double[atomn];
 res_be=new int*[resnum];
 for(i=0;i<resnum;i++)
	res_be[i]=new int[2];
 N=new double*[resnum];
 for(i=0;i<resnum;i++)
	N[i]=new double[3];
 CA=new double*[resnum];
 for(i=0;i<resnum;i++)
	CA[i]=new double[3];
 C=new double*[resnum];
 for(i=0;i<resnum;i++)
	C[i]=new double[3];
 CB=new double*[resnum];
 for(i=0;i<resnum;i++)
	CB[i]=new double[3];
 CB_miss=new int[resnum];

 for(i=0;i<atomn;i++)
	{
	 atomt[i]=m_c.atomt[i];
	 x[i]=m_c.x[i];
	 y[i]=m_c.y[i];
	 z[i]=m_c.z[i];
	}
 for(i=0;i<resnum;i++)
	{
	 res_t[i]=m_c.res_t[i];
	 res_chain[i]=m_c.res_chain[i];
	 res_be[i][0]=m_c.res_be[i][0];
	 res_be[i][1]=m_c.res_be[i][1];
	 res_id[i]=m_c.resid[m_c.res_be[i][0]];
	 N[i][0]=m_c.N[i][0];
	 N[i][1]=m_c.N[i][1];
	 N[i][2]=m_c.N[i][2];
	 CA[i][0]=m_c.CA[i][0];
	 CA[i][1]=m_c.CA[i][1];
	 CA[i][2]=m_c.CA[i][2];
	 C[i][0]=m_c.C[i][0];
	 C[i][1]=m_c.C[i][1];
	 C[i][2]=m_c.C[i][2];
	 CB[i][0]=m_c.CB[i][0];
	 CB[i][1]=m_c.CB[i][1];
	 CB[i][2]=m_c.CB[i][2];
	 CB_miss[i]=m_c.CB_miss[i];
	}
}

int Build_m::guess_CB_miss(Topprm *top)
{
 int i,j,k;
 double cd,bcd,abcd;
 Guess_d gd;

 for(i=0;i<resnum;i++)
 if(CB_miss[i]==1)
	{
	 for(j=0;j<top[i].icn;j++)
	 if(top[i].resic[j].atomid[3]=="CB")
		{
		 cd=top[i].resic[j].ic[4];
		 bcd=top[i].resic[j].ic[3];
		 abcd=top[i].resic[j].ic[2];
		 break;
		}
	 gd.init(N[i],C[i],CA[i],abcd,bcd,cd);
	 gd.guess();
	 gd.getd(CB[i]);
	}
}

int Build_m::build(Main_chain m_c,Topprm *top)
{
 get_mc(m_c);
 guess_CB_miss(top);
}

int COOR_m::init_cm(Build_m &b_m)
{
 int i,j;
 atomn=b_m.atomn;
 resnum=b_m.resnum;
 x=new double[atomn];
 y=new double[atomn];
 z=new double[atomn];
 N=new double*[resnum];
 CA=new double*[resnum];
 C=new double*[resnum];
 CB=new double*[resnum];
 for(i=0;i<resnum;i++)
	{
	 N[i]=new double[3];
	 CA[i]=new double[3];
	 C[i]=new double[3];
	 CB[i]=new double[3];
	}
 for(i=0;i<atomn;i++)
	{
	 x[i]=b_m.x[i];
	 y[i]=b_m.y[i];
	 z[i]=b_m.z[i];
	}
 for(i=0;i<resnum;i++)
 for(j=0;j<3;j++)
	{
	 N[i][j]=b_m.N[i][j];
	 CA[i][j]=b_m.CA[i][j];
	 C[i][j]=b_m.C[i][j];
	 CB[i][j]=b_m.CB[i][j];
	}
}

int COOR_s::operator=(COOR_s &c2)
{
 int i;
 atomn=c2.atomn;
 for(i=0;i<atomn;i++)
	{
	 x[i]=c2.x[i];
	 y[i]=c2.y[i];
	 z[i]=c2.z[i];
	}
 for(i=0;i<3;i++)
	{
	 N[i]=c2.N[i];
	 CA[i]=c2.CA[i];
	 C[i]=c2.C[i];
	 CB[i]=c2.CB[i];
	}
}

int Build_s::get_ss(Build_m m_c,Read_res r_r)
{
 int i,j,k,l;
 resnum=m_c.resnum;
 res_t=new string[resnum];
 ss=new Structure_s[resnum];

 for(i=0;i<resnum;i++)
	res_t[i]=m_c.res_t[i];

 for(i=0;i<resnum;i++)
 for(j=0;j<r_r.resn;j++)
 if(res_t[i]==r_r.res[j].res_t)
	{
	 ss[i].atomn=r_r.res[j].atomn_s;
	 ss[i].atomt=new string[ss[i].atomn];
	 for(k=0;k<ss[i].atomn;k++)
		ss[i].atomt[k]=r_r.res[j].atomt_s[k];
	 ss[i].axisn=r_r.res[j].axisn;
	 ss[i].axis=new int[ss[i].axisn];
	 ss[i].rttn=new int[ss[i].axisn];
	 ss[i].rtt=new int*[ss[i].axisn];
	 for(k=0;k<ss[i].axisn;k++)
		{
		 ss[i].axis[k]=r_r.res[j].axis[k];
		 ss[i].rttn[k]=r_r.res[j].rttn[k];
		 ss[i].rtt[k]=new int[ss[i].rttn[k]];
		 for(l=0;l<ss[i].rttn[k];l++)
			ss[i].rtt[k][l]=r_r.res[j].rtt[k][l];
		}
	 break;
	}
}

int Build_s::get_cs(Res_rtm rt)
{
 int i,j,k1,k2,k3;
 rtmn=new int[resnum];
 cs=new COOR_s*[resnum];
 for(i=0;i<resnum;i++)
	{
	 for(j=0;j<rt.rtn;j++)
	 if(res_t[i]==rt.rtmrt[j])
		{
		 rtmn[i]=rt.rtmn[j];
		 cs[i]=new COOR_s[rtmn[i]];
		 for(k1=0;k1<rtmn[i];k1++)
			{
			 cs[i][k1].x=new double[ss[i].atomn];
			 cs[i][k1].y=new double[ss[i].atomn];
			 cs[i][k1].z=new double[ss[i].atomn];
			 cs[i][k1].atomn=ss[i].atomn;
			 for(k2=0;k2<ss[i].atomn;k2++)
			 for(k3=0;k3<rt.rtmcr[j][k1].n;k3++)
			 if(ss[i].atomt[k2]==rt.rtmcr[j][k1].atomt[k3])
				{
				 cs[i][k1].x[k2]=rt.rtmcr[j][k1].x[k3];
				 cs[i][k1].y[k2]=rt.rtmcr[j][k1].y[k3];
				 cs[i][k1].z[k2]=rt.rtmcr[j][k1].z[k3];
				 break;
				}
			 cs[i][k1].N=new double[3];
			 cs[i][k1].CA=new double[3];
			 cs[i][k1].C=new double[3];
			 cs[i][k1].CB=new double[3];
			 for(k2=0;k2<3;k2++)
				cs[i][k1].N[k2]=rt.rtmcr[j][k1].N[k2];
			 for(k2=0;k2<3;k2++)
				cs[i][k1].CA[k2]=rt.rtmcr[j][k1].CA[k2];
			 for(k2=0;k2<3;k2++)
				cs[i][k1].C[k2]=rt.rtmcr[j][k1].C[k2];
			 for(k2=0;k2<3;k2++)
				cs[i][k1].CB[k2]=rt.rtmcr[j][k1].CB[k2];
			}
		 break;
		}
	}
}

int Build_s::rt_fit(Build_m m_c)
{
 int i,j,k;
 
 for(i=0;i<resnum;i++)
 for(j=0;j<rtmn[i];j++)
	{
	 if(cs[i][j].atomn>0)
	 fit_one(m_c,i,j);
	}
}

int Build_s::fit_one(Build_m m_c,int n1,int n2)
{
 int i,j,k;
 double theta;
 double m_v1,m_b1;
 double a1[3],b1[3],a2[3],b2[3],v1[3];
 double CA_t[3];

 for(i=0;i<3;i++)
	CA_t[i]=cs[n1][n2].CA[i];
 for(i=0;i<cs[n1][n2].atomn;i++)
	{
	 cs[n1][n2].x[i]-=CA_t[0];
	 cs[n1][n2].y[i]-=CA_t[1];
	 cs[n1][n2].z[i]-=CA_t[2];
	}
 for(i=0;i<3;i++)
	{
	 cs[n1][n2].N[i]-=CA_t[i];
	 cs[n1][n2].CA[i]-=CA_t[i];
	 cs[n1][n2].C[i]-=CA_t[i];
	 cs[n1][n2].CB[i]-=CA_t[i];
	}
 for(i=0;i<3;i++)
	{a1[i]=m_c.CB[n1][i]-m_c.CA[n1][i];
	 b1[i]=cs[n1][n2].CB[i]-cs[n1][n2].CA[i];
	}
 theta=get_angle(a1[0],a1[1],a1[2],b1[0],b1[1],b1[2]);
 L_CHAX(b1[0],b1[1],b1[2],a1[0],a1[1],a1[2],v1[0],v1[1],v1[2]);
 L_MO(v1[0],v1[1],v1[2],m_v1);
 v1[0]=v1[0]/m_v1; v1[1]=v1[1]/m_v1; v1[2]=v1[2]/m_v1;
 rotation_a(cs[n1][n2].x,cs[n1][n2].y,cs[n1][n2].z,cs[n1][n2].atomn,v1,theta);
 rotation(cs[n1][n2].N,v1,theta);
 rotation(cs[n1][n2].C,v1,theta);
 rotation(cs[n1][n2].CB,v1,theta);


 for(i=0;i<3;i++)
	{
	 b1[i]=cs[n1][n2].CB[i]-cs[n1][n2].CA[i];
	 a2[i]=m_c.CA[n1][i]-m_c.N[n1][i];
	 b2[i]=cs[n1][n2].N[i]-cs[n1][n2].CA[i];
	}
 theta=get_dhd(a2,b1,b2);
 L_MO(b1[0],b1[1],b1[2],m_b1);
 b1[0]=b1[0]/m_b1; b1[1]=b1[1]/m_b1; b1[2]=b1[2]/m_b1;

 rotation_a(cs[n1][n2].x,cs[n1][n2].y,cs[n1][n2].z,cs[n1][n2].atomn,b1,-theta);
 rotation(cs[n1][n2].N,b1,-theta);
 rotation(cs[n1][n2].C,b1,-theta);
 rotation(cs[n1][n2].CB,b1,-theta);

 for(i=0;i<cs[n1][n2].atomn;i++)
	{
	 cs[n1][n2].x[i]+=m_c.CA[n1][0];
	 cs[n1][n2].y[i]+=m_c.CA[n1][1];
	 cs[n1][n2].z[i]+=m_c.CA[n1][2];
	}
 for(i=0;i<3;i++)
	{
	 cs[n1][n2].N[i]+=m_c.CA[n1][i];
	 cs[n1][n2].CA[i]+=m_c.CA[n1][i];
	 cs[n1][n2].C[i]+=m_c.CA[n1][i];
	 cs[n1][n2].CB[i]+=m_c.CA[n1][i];
	}
}

int Build_s::build(Build_m m_c,Read_res r_r,Res_rtm rt)
{
 get_ss(m_c,r_r);
 get_cs(rt);
 rt_fit(m_c);
}

int Cpy_s::init_cs(Build_s &b_s)
{
 int i,n;
 n=b_s.resnum;
 cs=new COOR_s[n];
 for(i=0;i<n;i++)
	{
	 cs[i].atomn=b_s.ss[i].atomn;
	 cs[i].x=new double[cs[i].atomn];
	 cs[i].y=new double[cs[i].atomn];
	 cs[i].z=new double[cs[i].atomn];
	 cs[i].N=new double[3];
	 cs[i].CA=new double[3];
	 cs[i].C=new double[3];
	 cs[i].CB=new double[3];
	}
}
