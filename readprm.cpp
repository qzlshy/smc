#include "long.h"

int Read_prm::get_prm(ifstream *f)
{
 int i;
 string l;
 string st;
 stringstream s;

 f->clear();
 f->seekg(0);
 while(!f->eof())
	{
	 getline(*f,l);
	 s.clear();
         s.str(l);
	 s>>st;
	 if(st=="BONDS") break;
	}
 nbond=0;
 while(!f->eof())
	{getline(*f,l);
	 if(l.length()==0) continue;
	 s.clear();
         s.str(l);
	 s>>st;
	 if(st=="ANGLES"||st=="END") break;
	 if(st[0]=='!'||st=="\t") continue;
	 nbond++;
	}

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="ANGLES") break;
        }
 nangle=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="DIHEDRALS"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         nangle++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="DIHEDRALS") break;
        }
 ndihedral=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="IMPROPER"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         ndihedral++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        { 
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="IMPROPER") break;
        }
 nimpro=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="CMAP"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         nimpro++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="NONBONDED") break;
        }
 getline(*f,l);
 nnonbond=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="HBOND"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         nnonbond++;
        }

 bond=new BONDS[nbond];
 angle=new ANGLES[nangle];
 dihedral=new DIHEDRALS[ndihedral];
 improp=new IMPROPER[nimpro];
 nonbond=new NONBONDED[nnonbond];

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="BONDS") break;
        }

 i=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="ANGLES"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
	 bond[i].a1=st;
	 s>>bond[i].a2>>bond[i].Kb>>bond[i].b0;
	 i++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="ANGLES") break;
        }

 i=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="DIHEDRALS"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         angle[i].a1=st;
         s>>angle[i].a2>>angle[i].a3>>angle[i].Ktheta>>angle[i].Theta0>>angle[i].Kub>>angle[i].S0;
         i++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="DIHEDRALS") break;
        }

 i=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="IMPROPER"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         dihedral[i].a1=st;
         s>>dihedral[i].a2>>dihedral[i].a3>>dihedral[i].a4>>dihedral[i].Kchi>>dihedral[i].n>>dihedral[i].delta;
         i++;
        }


 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="IMPROPER") break;
        }

 i=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="CMAP"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         improp[i].a1=st;
         s>>improp[i].a2>>improp[i].a3>>improp[i].a4>>improp[i].Kpsi>>improp[i].n>>improp[i].psi0;
         i++;
        }

 f->clear();
 f->seekg(0);
 while(!f->eof())
        {
         getline(*f,l);
         s.clear();
         s.str(l);
         s>>st;
         if(st=="NONBONDED") break;
        }
 getline(*f,l);
 i=0;
 while(!f->eof())
        {getline(*f,l);
         if(l.length()==0) continue;
         s.clear();
         s.str(l);
         s>>st;
         if(st=="HBOND"||st=="END") break;
         if(st[0]=='!'||st=="\t") continue;
         nonbond[i].a1=st;
         s>>nonbond[i].ignored1>>nonbond[i].epsilon>>nonbond[i].Rmin2>>nonbond[i].ignored2>>nonbond[i].eps14>>nonbond[i].Rmin214;
         i++;
        }

}
