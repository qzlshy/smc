#include "long.h"

int Read_res::getres(ifstream *f)
{
 int i,j;
 string tmp;
 resn=22;
 res=new RES_SIDE[22];
 for(i=0;i<22;i++)
	{
	 *f>>tmp>>res[i].res_t;
	 *f>>tmp>>res[i].atomn;
	 res[i].atomt=new string[res[i].atomn];
	 for(j=0;j<res[i].atomn;j++)
		*f>>res[i].atomt[j];
	 *f>>tmp>>res[i].axisn;
	 res[i].rotate=new string*[res[i].axisn];
	 for(j=0;j<res[i].axisn;j++)
		res[i].rotate[j]=new string[4];
	 for(j=0;j<res[i].axisn;j++)
		*f>>res[i].rotate[j][0]>>res[i].rotate[j][1]>>res[i].rotate[j][2]>>res[i].rotate[j][3];
	}
}

int Read_res::get_nh()
{
 int i,j,k;
 for(i=0;i<22;i++)
	{
	 k=0;
	 for(j=0;j<res[i].atomn;j++)
		{
		 if(res[i].atomt[j][0]!='H')
			{
			 res[i].atomt[k]=res[i].atomt[j];
			 k++;
			}
		} 
	 res[i].atomn=k;
	}
}

int Read_res::getside()
{
 int i,j,k,k2,t;
 int n1,n2;

 string mc_atomt[15]={"N","HN","CA","HA","C","O","HA1","HA2","HT1","HT2","HT3","HN1","HN2","OT1","OT2"};
 string bgdezh="BGDEZH";

 for(i=0;i<22;i++)
	{
	 res[i].atomn_s=0;
	 for(j=0;j<res[i].atomn;j++)
		{
		 t=0;
		 for(k=0;k<15;k++)
		 if(res[i].atomt[j]==mc_atomt[k])
			{
			 t=1;
			 break;
			}
		 if(t==0)
			res[i].atomn_s++;
		}
	 res[i].atomt_s=new string[res[i].atomn_s];
	}

 for(i=0;i<22;i++)
        {
         res[i].atomn_s=0;
         for(j=0;j<res[i].atomn;j++)
                {
                 t=0;
                 for(k=0;k<15;k++)
                 if(res[i].atomt[j]==mc_atomt[k])
                        {
                         t=1;
                         break;
                        }
                 if(t==0)
                        {
			 res[i].atomt_s[res[i].atomn_s]=res[i].atomt[j];
			 res[i].atomn_s++;
			}
                }
        }

 for(i=0;i<22;i++)
	{
	 res[i].axis=new int[res[i].axisn];
	 res[i].rttn=new int[res[i].axisn];
	 res[i].rtt=new int*[res[i].axisn];
	 for(j=0;j<res[i].axisn;j++)
		res[i].rtt[j]=new int[res[i].atomn_s];
	 for(j=0;j<res[i].axisn;j++)
	 	{
		 for(k=0;k<res[i].atomn_s;k++)
		 if(res[i].atomt_s[k]==res[i].rotate[j][2])
			{res[i].axis[j]=k;
			 break;
			}
		}
	 for(j=0;j<res[i].axisn;j++)
		{
		 for(k=0;k<6;k++)
		 if(res[i].rotate[j][2][1]==bgdezh[k])
			{
			 n1=k; break;
			}
		 t=0;
		 if(res[i].rotate[j][2].size()==2)
			{for(k=0;k<res[i].atomn_s;k++)
				{
				 for(k2=0;k2<6;k2++)
				 if(res[i].atomt_s[k][1]==bgdezh[k2])
					{n2=k2; break;
					}
				 if(n2>=n1)
					{
					 res[i].rtt[j][t]=k;
					 t++;
					}
				}
			 res[i].rttn[j]=t;
			}
		 t=0;
		 if(res[i].rotate[j][2].size()==3)
			{for(k=0;k<res[i].atomn_s;k++)
                                {
                                 for(k2=0;k2<6;k2++)
                                 if(res[i].atomt_s[k][1]==bgdezh[k2])
                                        {n2=k2; break;
                                        }
                                 if(n2>n1||(n2==n1&&res[i].atomt_s[k][2]==res[i].rotate[j][2][2]))
                                        {
                                         res[i].rtt[j][t]=k;
                                         t++;
                                        }
                                }
                         res[i].rttn[j]=t;
                        }
		}
	}
}

int Read_res::get_res(ifstream *fp1)
{
 getres(fp1);
 get_nh();
 getside();
}


int Read_res::get_res(ifstream *fp1,int useh)
{
 getres(fp1);
 if(useh==0) get_nh();
 getside();
}
