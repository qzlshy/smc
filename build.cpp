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

int Build_s::get_ss(Main_chain m_c,Read_res r_r)
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
