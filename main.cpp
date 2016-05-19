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
}

