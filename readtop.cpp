#include "long.h"

int Topprm::getan(ifstream *infile)
{
 string tmps;
 int i;
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="GROUP")
                infile->ignore(1024,'\n');
         else if(tmps=="ATOM")
                {i++;
                 infile->ignore(1024,'\n');
                }
         else break;
        }
 return(i);
}

int Topprm::getin(ifstream *infile)
{
 string tmps;
 int i;
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps!="IC")
                infile->ignore(1024,'\n');
         else
                {i++;
                 infile->ignore(1024,'\n');
                 break;
                }
        }

 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC")
                {infile->ignore(1024,'\n');
                 i++;
                }
         else break;
        }
 return(i);
}

int Topprm::readatom(ifstream *infile)
{
 string tmps;
 int i;

 atom=new string[atomn];
 atom2=new string[atomn];
 charge=new double[atomn];

 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="GROUP")
                infile->ignore(1024,'\n');
         else if(tmps=="ATOM")
                {(*infile)>>atom[i]>>atom2[i]>>charge[i];
                 infile->ignore(1024,'\n');
                 i++;
                }
         else break;
        }
 return(i);
}

int Topprm::readic(ifstream *infile)
{
 string tmps;
 int i;

 resic=new Topic[icn];

 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps!="IC")
                infile->ignore(1024,'\n');
         else
                {(*infile)>>resic[i].atomid[0]>>resic[i].atomid[1]>>resic[i].atomid[2]>>resic[i].atomid[3]>>resic[i].ic[0]>>resic[i].ic[1]>>resic[i].ic[2]>>resic[i].ic[3]>>resic[i].ic[4];
                 infile->ignore(1024,'\n');
                 i++;
                 break;
                }
        }

 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC")
                {(*infile)>>resic[i].atomid[0]>>resic[i].atomid[1]>>resic[i].atomid[2]>>resic[i].atomid[3]>>resic[i].ic[0]>>resic[i].ic[1]>>resic[i].ic[2]>>resic[i].ic[3]>>resic[i].ic[4];
                 infile->ignore(1024,'\n');
                 i++;
                }
         else break;
        }
 return(i);
}

int Topprm::getimn(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC") break;
         if(tmps=="IMPR")
                {
                 while(1)
                        {ch=infile->peek();
                         if(ch==' ')
                                ch=infile->get();
                         else if(ch=='\n')
                                break;
                         else
                                {
                                 (*infile)>>tmps>>tmps>>tmps>>tmps;
                                 i++;
                                }
                        }
                }
        }
 return(i);
}

int Topprm::readimpr(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 impr=new string*[imprn];
 for(i=0;i<imprn;i++)
	impr[i]=new string[4];
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC") break;
         if(tmps=="IMPR")
                {
                 while(1)
                        {ch=infile->peek();
                         if(ch==' ')
                                ch=infile->get();
                         else if(ch=='\n')
                                break;
                         else
                                {
                                 (*infile)>>impr[i][0]>>impr[i][1]>>impr[i][2]>>impr[i][3];
                                 i++;
                                }
                        }
                }
        }
 return(i);
}

int Topprm::getcn(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC") break;
         if(tmps=="CMAP")
                {
                 while(1)
                        {ch=infile->peek();
                         if(ch==' ')
                                ch=infile->get();
                         else if(ch=='\n')
                                break;
                         else
                                {
                                 (*infile)>>tmps>>tmps>>tmps>>tmps;
                                 i++;
                                }
                        }
                }
        }
 return(i);
}

int Topprm::readcmap(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 CMAP=new string*[cmapn];
 for(i=0;i<cmapn;i++)
	CMAP[i]=new string[4];
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC") break;
         if(tmps=="CMAP")
                {
                 while(1)
                        {ch=infile->peek();
                         if(ch==' ')
                                ch=infile->get();
                         else if(ch=='\n')
                                break;
                         else
                                {
                                 (*infile)>>CMAP[i][0]>>CMAP[i][1]>>CMAP[i][2]>>CMAP[i][3];
                                 i++;
                                }
                        }
                }
        }
 return(i);
}

int Topprm::getbn(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 i=0;
 while(1)
	{(*infile)>>tmps;
	 if(tmps=="IC") break;
	 if(tmps=="BOND"||tmps=="DOUBLE")
		{
		 while(1)
			{ch=infile->peek();
			 if(ch==' ')
				ch=infile->get();
			 else if(ch=='\n')
				break;
			 else
				{
				 (*infile)>>tmps;
				 i++;
				}
			}
		}
	}
 return(i/2);
}

int Topprm::readbond(ifstream *infile)
{
 string tmps;
 int i,k;
 char ch;
 bond=new string*[bondn];
 for(i=0;i<bondn;i++)
	bond[i]=new string[2];
 i=0;
 while(1)
        {(*infile)>>tmps;
         if(tmps=="IC") break;
         if(tmps=="BOND"||tmps=="DOUBLE")
                {
                 while(1)
                        {ch=infile->peek();
                         if(ch==' ')
                                ch=infile->get();
                         else if(ch=='\n')
                                break;
                         else
                                {
                                 (*infile)>>bond[i][0]>>bond[i][1];
                                 i++;
                                }
                        }
                }
        }
}

int Topprm::readtop(ifstream *infile,string res)
{
 int i;
 string tmps;
 streampos pos;

 infile->clear();
 infile->seekg(0);
 while(1)
        {(*infile)>>tmps;
         if(tmps=="RESI")
                {(*infile)>>tmps;
                 if(tmps==res)
                        {resid=tmps;
                         infile->ignore(1024,'\n');
                         break;
                        }
                }
        }
 pos=(*infile).tellg();
 atomn=getan(infile);
 infile->clear();
 infile->seekg(pos);
 icn=getin(infile);
 infile->clear();
 infile->seekg(pos);
 bondn=getbn(infile);
 infile->clear();
 infile->seekg(pos);
 imprn=getimn(infile);
 infile->clear();
 infile->seekg(pos);
 cmapn=getcn(infile);
 infile->clear();
 infile->seekg(pos);

 readatom(infile);
 infile->clear();
 infile->seekg(pos);
 readic(infile);
 infile->clear();
 infile->seekg(pos);
 readbond(infile);
 infile->clear();
 infile->seekg(pos);
 readimpr(infile);
 infile->clear();
 infile->seekg(pos);
 readcmap(infile);
 infile->clear();
 infile->seekg(0);


 return(1);
}

int Topprm::readtopter(ifstream *infile,string res)
{
 int i;
 string tmps;
 streampos pos;

 infile->clear();
 infile->seekg(0);
 while(1)
        {(*infile)>>tmps;
         if(tmps=="PRES")
                {(*infile)>>tmps;
                 if(tmps==res)
                        {resid=tmps;
                         infile->ignore(1024,'\n');
                         break;
                        }
                }
        }
 pos=(*infile).tellg();
 atomn=getan(infile);
 infile->clear();
 infile->seekg(pos);
 icn=getin(infile);
 infile->clear();
 infile->seekg(pos);
 bondn=getbn(infile);
 infile->clear();
 infile->seekg(pos);

 readatom(infile);
 infile->clear();
 infile->seekg(pos);
 readic(infile);
 infile->clear();
 infile->seekg(pos);
 readbond(infile);
 infile->clear();
 infile->seekg(0);


 return(1);
}
