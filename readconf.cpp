#include "long.h"

Read_conf::Read_conf(char *s)
{
 int i,j,k;
 ifstream fp(s);
 string t;
 k=0;
 while(1)
        {
         getline(fp,t);
         if(fp.eof()) break;
         k++;
        }
 n=k;
 line=new string[n];
 fp.clear();
 fp.seekg(0);

 for(i=0;i<n;i++)
        {
         getline(fp,line[i]);
        }
 fp.close();
}

int Read_conf::inputconf(char *s)
{
 int i,j,k;
 ifstream fp(s);
 string t;
 k=0;
 while(1)
	{
	 getline(fp,t);
	 if(fp.eof()) break;
	 k++;
	}
 n=k;
 line=new string[n];
 fp.clear();
 fp.seekg(0);

 for(i=0;i<n;i++)
	{
	 getline(fp,line[i]);
	}
 fp.close();
}

string Read_conf::getstr(char *s)
{
 int i,j,k;
 string t1,t2;
 stringstream ss;

 for(i=0;i<n;i++)
	{
	 ss.str(line[i]);
	 ss>>t1;
	 if(t1==s)
		{
		 ss>>t2;
		 ss.clear();
		 return(t2);
		}
	 ss.clear();
	}
 cout<<"No "<<s<<" exist in conf file!!!\n";
 exit(0);
}

int Read_conf::getint(char *s)
{
 int i,j,k;
 string t1,t2;
 stringstream ss;

 for(i=0;i<n;i++)
        {
         ss.str(line[i]);
         ss>>t1;
         if(t1==s)
                {
                 ss>>k;
                 ss.clear();
                 return(k);
                }
         ss.clear();
        }
 cout<<"No "<<s<<" exist in conf file!!!\n";
 exit(0);
}

int Read_conf::getfloat(char *s)
{
 int i,j;
 float k;
 string t1,t2;
 stringstream ss;

 for(i=0;i<n;i++)
        {
         ss.str(line[i]);
         ss>>t1;
         if(t1==s)
                {
                 ss>>k;
                 ss.clear();
                 return(k);
                }
         ss.clear();
        }
 cout<<"No "<<s<<" exist in conf file!!!\n";
 exit(0);
}
