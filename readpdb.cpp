#include "long.h"

int Read_pdb::readpdb(ifstream *pdbfile)
{
 string l,l1,tmps,tmpres;
 stringstream stream;
 int i,k,n;
 n=0;
 while(1)
	{getline((*pdbfile),l);
	 if(pdbfile->eof()) break;
	 if(l.empty()) continue;
	 if(l.size()<4) continue;
	 if(l.compare(0,4,"ATOM")) continue;
	 stream.clear();
	 stream.str(l);
	 stream>>l1;
	 if(l1=="ATOM")
		n++;
	}
 atomn=n;
 pdbfile->clear();
 pdbfile->seekg(0);
 atomt=new string[atomn];
 alt=new char[atomn];
 residt=new string[atomn];
 chain=new char[atomn];
 serial=new int[atomn];
 resid= new int[atomn];
 index=new int[atomn];
 segid= new string[atomn];
 x=new double[atomn];
 y=new double[atomn];
 z=new double[atomn];
 tof=new double[atomn];
 line=new string[atomn];

 n=0;
 while(1)
        {
         getline((*pdbfile),l);
         if(pdbfile->eof()) break;
         if(l.empty()) continue;
         if(l.size()<4) continue;
         if(l.compare(0,4,"ATOM")) continue;
	 line[n]=l;
         stream.clear();
         stream.str(l);
         stream>>l1;
         if(l1=="ATOM") {
         tmps=l.substr(6,5);
         stream.clear();
         stream.str(tmps);
         stream>>serial[n];
         tmps=l.substr(12,4);
         stream.clear();
         stream.str(tmps);
         stream>>atomt[n];
         tmps=l.substr(16,1);
         stream.clear();
         stream.str(tmps);
         stream>>alt[n];
         tmps=l.substr(17,3);
         stream.clear();
         stream.str(tmps);
         stream>>residt[n];
         tmps=l.substr(21,1);
         stream.clear();
         stream.str(tmps);
         stream>>chain[n];
         tmps=l.substr(22,4);
         stream.clear();
         stream.str(tmps);
         stream>>resid[n];
         tmps=l.substr(30,8);
         stream.clear();
         stream.str(tmps);
         stream>>x[n];
         tmps=l.substr(38,8);
         stream.clear();
         stream.str(tmps);
         stream>>y[n];
         tmps=l.substr(46,8);
         stream.clear();
         stream.str(tmps);
         stream>>z[n];
         tmps=l.substr(60,6);
         stream.clear();
         stream.str(tmps);
         stream>>tof[n];
	 if(l.size()<76) continue;
	 tmps=l.substr(72,4);
	 stream.clear();
	 stream.str(tmps);
	 stream>>segid[n];
         index[n]=n;
         n++; }
        }

 

 pdbfile->clear();
 pdbfile->seekg(0);
}

int Read_pdb::getprotein()
{
 int i,j,k,n;
 string acids[25]= {"ALA","ARG","ASN","ASP","CYS", \
			"GLN","GLU","GLY","HSD","HSE", \
			"HSP","ILE","LEU","LYS","MET", \
			"PHE","PRO","SER","THR","TRP", \
			"TYR","VAL","HIS","HIE","MSE"};
 int acidsn=25;
 n=0;
 for(i=0;i<atomn;i++)
	{k=0;
	 for(j=0;j<acidsn;j++)
		{
		 if(residt[i]==acids[j])
			{
			 k=1;
			 break;
			}
		}
	 if(k==1)
		{
		 serial[n]=serial[i];
		 resid[n]=resid[i];
		 index[n]=index[i];
		 atomt[n]=atomt[i];
		 residt[n]=residt[i];
		 segid[n]=segid[i];
		 tof[n]=tof[i];
		 alt[n]=alt[i];
		 chain[n]=chain[i];
		 x[n]=x[i];
		 y[n]=y[i];
		 z[n]=z[i];
		 line[n]=line[i];
		 n++;
		}
	}
 atomn=n;
}

int Read_pdb::choosealt()
{
 int i,n;
 char talt;

 for(i=0;i<atomn;i++)
 if(alt[i]!='\0')
	{
	 talt=alt[i];
	 break;
	}
 n=0;
 for(i=0;i<atomn;i++)
 if(alt[i]=='\0'||alt[i]==talt)
	{
	 serial[n]=serial[i];
	 resid[n]=resid[i];
	 index[n]=index[i];
	 atomt[n]=atomt[i];
	 residt[n]=residt[i];
	 segid[n]=segid[i];
	 tof[n]=tof[i];
	 alt[n]=alt[i];
	 chain[n]=chain[i];
	 x[n]=x[i];
	 y[n]=y[i];
	 z[n]=z[i];
	 line[n]=line[i];
	 n++;
	}
 atomn=n;
}

int Read_pdb::removeh()
{
 int i,n;

 n=0;
 for(i=0;i<atomn;i++)
 if(atomt[i][0]!='H')
	{
	 serial[n]=serial[i];
	 resid[n]=resid[i];
	 index[n]=index[i];
	 atomt[n]=atomt[i];
	 residt[n]=residt[i];
	 segid[n]=segid[i];
	 tof[n]=tof[i];
	 alt[n]=alt[i];
	 chain[n]=chain[i];
	 x[n]=x[i];
	 y[n]=y[i];
	 z[n]=z[i];
	 line[n]=line[i];
	 n++;
	}
 atomn=n;
}

int Pdb_res::getres()
{
 int i,k,n;
 n=1;
 for(i=1;i<atomn;i++)
 if(resid[i]!=resid[i-1]||residt[i]!=residt[i-1]||chain[i]!=chain[i-1])
	n++;
 resnum=n;
 res_be=new int*[resnum];
 for(i=0;i<resnum;i++)
	res_be[i]=new int[2];

 n=0;
 res_be[n][0]=0;
 for(i=1;i<atomn;i++)
 if(resid[i]!=resid[i-1]||residt[i]!=residt[i-1]||chain[i]!=chain[i-1])
	{
	 res_be[n][1]=i;
	 n++;
	 res_be[n][0]=i;
	}
 res_be[n][1]=i;

 res_t=new string[resnum];
 res_chain=new char[resnum];
 for(i=0;i<resnum;i++)
	{res_t[i]=residt[res_be[i][0]];
	 res_chain[i]=chain[res_be[i][0]];
	}

 return(resnum); 
}

int Res_rtm::get_rtm(ifstream *file)
{
 int i,j,k,l,m,n;
 string rttmp;


 readpdb(file);
 getres();

 rttmp=res_t[0]; k=1;
 for(i=0;i<resnum;i++)
 if(rttmp!=res_t[i])
	{
	 rttmp=res_t[i];
	 k++;
	}
 rtn=k;
 rtmn=new int[rtn];
 rtmrt=new string[rtn];

 k=0;
 rtmrt[k]=res_t[0];
 rttmp=res_t[0];
 k++;
 for(i=0;i<resnum;i++)
 if(rttmp!=res_t[i])
	{
	 rtmrt[k]=res_t[i];
	 rttmp=res_t[i];
	 k++;
	}
 rtmcr=new COOR*[rtn];
 for(i=0;i<rtn;i++)
	{
	 k=0;
	 for(j=0;j<resnum;j++)
	 if(res_t[j]==rtmrt[i])
		{
		 k++;
		}
	 rtmn[i]=k;
	 rtmcr[i]=new COOR[rtmn[i]];
	 k=0;
	 for(j=0;j<resnum;j++)
	 if(res_t[j]==rtmrt[i])
		{
		 m=res_be[j][1]-res_be[j][0];
		 rtmcr[i][k].n=m;
		 rtmcr[i][k].atomt=new string[m];
		 rtmcr[i][k].x=new double[m];
		 rtmcr[i][k].y=new double[m];
		 rtmcr[i][k].z=new double[m];
		 n=0;
		 for(l=res_be[j][0];l<res_be[j][1];l++)
			{
			 rtmcr[i][k].atomt[n]=atomt[l];
			 rtmcr[i][k].x[n]=x[l];
			 rtmcr[i][k].y[n]=y[l];
			 rtmcr[i][k].z[n]=z[l];
			 n++;
			}
		 for(l=0;l<m;l++)
			{
			 if(rtmcr[i][k].atomt[l]=="N")
				{
				 rtmcr[i][k].N[0]=rtmcr[i][k].x[l];
				 rtmcr[i][k].N[1]=rtmcr[i][k].y[l];
				 rtmcr[i][k].N[2]=rtmcr[i][k].z[l];
				}
			 if(rtmcr[i][k].atomt[l]=="CA")
				{
				 rtmcr[i][k].CA[0]=rtmcr[i][k].x[l];
				 rtmcr[i][k].CA[1]=rtmcr[i][k].y[l];
				 rtmcr[i][k].CA[2]=rtmcr[i][k].z[l];
				}
			 if(rtmcr[i][k].atomt[l]=="C")
				{
				 rtmcr[i][k].C[0]=rtmcr[i][k].x[l];
				 rtmcr[i][k].C[1]=rtmcr[i][k].y[l];
				 rtmcr[i][k].C[2]=rtmcr[i][k].z[l];
				}
			 if(rtmcr[i][k].atomt[l]=="CB"||rtmcr[i][k].atomt[l]=="HA1")
				{
				 rtmcr[i][k].CB[0]=rtmcr[i][k].x[l];
				 rtmcr[i][k].CB[1]=rtmcr[i][k].y[l];
				 rtmcr[i][k].CB[2]=rtmcr[i][k].z[l];
				}
			}
		 k++;
		}

	}
}

int Main_chain::get_mc(ifstream *fp1)
{
 int mc_n=15;
 string mc_atomt[15]={"N","HN","CA","HA","C","O","HA1","HA2","HT1","HT2","HT3","HN1","HN2","OT1","OT2"};
 int i,j,k,n;

 readpdb(fp1);
 getprotein();
 choosealt();
 removeh();
 getres();

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

 N_n=new int[resnum];
 CA_n=new int[resnum];
 C_n=new int[resnum];
 CB_miss=new int[resnum];

 for(i=0;i<resnum;i++)
        CB_miss[i]=1;

 for(i=0;i<resnum;i++)
 for(j=res_be[i][0];j<res_be[i][1];j++)
	{
	 if(atomt[j]=="N")
		{
		 N[i][0]=x[j];
		 N[i][1]=y[j];
		 N[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="CA")
		{
		 CA[i][0]=x[j];
		 CA[i][1]=y[j];
		 CA[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="C")
		{
		 C[i][0]=x[j];
		 C[i][1]=y[j];
		 C[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="CB"||atomt[j]=="HA1")
		{
		 CB[i][0]=x[j];
		 CB[i][1]=y[j];
		 CB[i][2]=z[j];
		 CB_miss[i]=0;
		 continue;
		}
	}

 n=0;
 for(i=0;i<atomn;i++)
	{
	 for(j=0;j<mc_n;j++)
	 if(mc_atomt[j]==atomt[i])
		{
		 serial[n]=serial[i];
		 resid[n]=resid[i];
		 index[n]=index[i];
		 atomt[n]=atomt[i];
		 residt[n]=residt[i];
		 segid[n]=segid[i];
		 tof[n]=tof[i];
		 alt[n]=alt[i];
		 chain[n]=chain[i];
		 x[n]=x[i];
		 y[n]=y[i];
		 z[n]=z[i];
		 line[n]=line[i];
		 n++;
		 break;
		}
	}
 atomn=n;

 n=0;
 res_be[n][0]=0;
 for(i=1;i<atomn;i++)
 if(resid[i]!=resid[i-1]||residt[i]!=residt[i-1]||chain[i]!=chain[i-1])
	{
	 res_be[n][1]=i;
	 n++;
	 res_be[n][0]=i;
	}
 res_be[n][1]=i;

 for(i=0;i<resnum;i++)
	{k=0;
	 for(j=res_be[i][0];j<res_be[i][1];j++)
		{
		 if(atomt[j]=="N")
			{
			 N_n[i]=j;
			 k++;
			 continue;
			}
		 if(atomt[j]=="CA")
			{
			 CA_n[i]=j;
			 k++;
			 continue;
			}
		 if(atomt[j]=="C")
			{
			 C_n[i]=j;
			 k++;
			 continue;
			}
		}
	 if(k!=3)
		{
		 cout<<"Main chain atom missing!!!\n";
		 exit(0);
		}
	}

 for(i=0;i<resnum;i++)
	{
	 if(res_t[i]=="HIS")
		res_t[i]="HSD";
	}

 return(resnum);
}


int Main_chain::get_mc(ifstream *fp1,int useh)
{
 int mc_n=15;
 string mc_atomt[15]={"N","HN","CA","HA","C","O","HA1","HA2","HT1","HT2","HT3","HN1","HN2","OT1","OT2"};
 int i,j,k,n;

 readpdb(fp1);
 getprotein();
 choosealt();
 if(useh==0) removeh();
 getres();

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

 N_n=new int[resnum];
 CA_n=new int[resnum];
 C_n=new int[resnum];
 CB_miss=new int[resnum];

 for(i=0;i<resnum;i++)
        CB_miss[i]=1;

 for(i=0;i<resnum;i++)
 for(j=res_be[i][0];j<res_be[i][1];j++)
	{
	 if(atomt[j]=="N")
		{
		 N[i][0]=x[j];
		 N[i][1]=y[j];
		 N[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="CA")
		{
		 CA[i][0]=x[j];
		 CA[i][1]=y[j];
		 CA[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="C")
		{
		 C[i][0]=x[j];
		 C[i][1]=y[j];
		 C[i][2]=z[j];
		 continue;
		}
	 if(atomt[j]=="CB"||atomt[j]=="HA1")
		{
		 CB[i][0]=x[j];
		 CB[i][1]=y[j];
		 CB[i][2]=z[j];
		 CB_miss[i]=0;
		 continue;
		}
	}

 n=0;
 for(i=0;i<atomn;i++)
	{
	 for(j=0;j<mc_n;j++)
	 if(mc_atomt[j]==atomt[i])
		{
		 serial[n]=serial[i];
		 resid[n]=resid[i];
		 index[n]=index[i];
		 atomt[n]=atomt[i];
		 residt[n]=residt[i];
		 segid[n]=segid[i];
		 tof[n]=tof[i];
		 alt[n]=alt[i];
		 chain[n]=chain[i];
		 x[n]=x[i];
		 y[n]=y[i];
		 z[n]=z[i];
		 line[n]=line[i];
		 n++;
		 break;
		}
	}
 atomn=n;

 n=0;
 res_be[n][0]=0;
 for(i=1;i<atomn;i++)
 if(resid[i]!=resid[i-1]||residt[i]!=residt[i-1]||chain[i]!=chain[i-1])
	{
	 res_be[n][1]=i;
	 n++;
	 res_be[n][0]=i;
	}
 res_be[n][1]=i;

 for(i=0;i<resnum;i++)
	{k=0;
	 for(j=res_be[i][0];j<res_be[i][1];j++)
		{
		 if(atomt[j]=="N")
			{
			 N_n[i]=j;
			 k++;
			 continue;
			}
		 if(atomt[j]=="CA")
			{
			 CA_n[i]=j;
			 k++;
			 continue;
			}
		 if(atomt[j]=="C")
			{
			 C_n[i]=j;
			 k++;
			 continue;
			}
		}
	 if(k!=3)
		{
		 cout<<"Main chain atom missing!!!\n";
		 exit(0);
		}
	}
 for(i=0;i<resnum;i++)
	{
	 if(res_t[i]=="HIS")
		res_t[i]="HSD";
	}

 return(resnum);
}
