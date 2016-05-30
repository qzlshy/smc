#include "long.h"

int main(int argc,char **argv)
{
 int i,j,k,resnum,m;
 Res_rtm r_t;
 Main_chain pdb_m;
 Read_res res;
 Topprm *top,nter,cter;
 Read_prm prm;
 
 Read_conf conf(argv[1]);
 ifstream ifrt(conf.getstr("rotamer").c_str());
 ifstream ifmc(conf.getstr("pdb").c_str());
 ifstream ifres(conf.getstr("res").c_str());
 ifstream iftop(conf.getstr("top").c_str());
 ifstream ifprm(conf.getstr("prm").c_str());
 m=conf.getint("copy");

 r_t.get_rtm(&ifrt);
 resnum=pdb_m.get_mc(&ifmc);
 res.get_res(&ifres);
 top=new Topprm[resnum];
 for(i=0;i<resnum;i++)
	top[i].readtop(&iftop,pdb_m.res_t[i]);
 nter.readtopter(&iftop,"NTER");
 cter.readtopter(&iftop,"CTER");
 prm.get_prm(&ifprm);

 Build_m b_m;
 b_m.build(pdb_m,top);

 Build_s b_s;
 b_s.build(b_m,res,r_t);

 COOR_m *c_m;
 c_m=new COOR_m[m];
 for(i=0;i<m;i++)
	c_m[i].init_cm(b_m);

 Cpy_s *cy_s;
 cy_s=new Cpy_s[m];
 for(i=0;i<m;i++)
        cy_s[i].init_cs(b_s);

 Rtm_tmp rt_t;
 rt_t.init(b_s);

 Caldist cdt;
 cdt.init(b_m,b_s);
 cdt.getw(0.7);

//now calculate



 srand(unsigned(time(0)));
 int rdi,t1,t2,t3;
 int rsn,*rs,rmn,*rm;
 double rd,*w,wt,result;
 w=new double[m];
 for(i=0;i<m;i++)
	w[i]=1.0;
 rs=new int[m];
 rm=new int[m];

 for(i=0;i<resnum;i++)
	{rsn=0; rmn=0;
	 for(j=0;j<m;j++)
		{
		 cdt.cal(c_m[j],cy_s[j],rt_t,i);
		 t1=cdt.goodn[i];
		 w[j]*=double(t1);
		 if(w[j]==0)
			{
			 rm[rmn]=j;
			 rmn++;
			}
		 else
			{
			 rdi=rand()%t1;
			 t2=cdt.good[i][rdi];
			 cy_s[j].cs[i]=rt_t.cs[i][t2];
			 cy_s[j].ipk[i]=1;
			 rs[rsn]=j;
			 rsn++;
			}
		}
	 if(rsn==0)
		{
		 cout<<"Not well!!!\n";
		 exit(0);
		}
	 for(j=0;j<rmn;j++)
		{
		 t1=rand()%rsn;
		 t2=rs[t1];
		 t3=rm[j];
		 cy_s[t3]=cy_s[t2];
		 w[t2]/=2.0;
		 w[t3]=w[t2];
		}
	}
 wt=0.0;
 for(i=0;i<m;i++)
	wt+=w[i];
 result=log(wt/double(m));
 cout<<result<<'\n';

 
 
}

