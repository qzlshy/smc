# smc
Calculate protein side chain conformational entropy for a given protein backbone conformation.<br/>

Compile<br/>
make clean<br/>
make<br/>

Run<br/>
./long cal.conf<br/>
<br/>
cal.conf is the configuration file, defind the input files of the program. <br/>

rotamer ./usefile/rotamer.pdb &nbsp;&nbsp;&nbsp;&nbsp; #rotamer library   <br/>
pdb ./ubq/ubq-1.pdb          &nbsp;&nbsp;&nbsp;&nbsp; #input protein structure, you can change to yourself file<br/>
res ./usefile/res.txt           &nbsp;&nbsp;&nbsp;&nbsp;           #defind side chain dihedral angle of 20 type residues  <br/>
top ./usefile/top_all36_prot_lipid.rtf   &nbsp;&nbsp;&nbsp;&nbsp;;  #Topology file  <br/>
prm ./usefile/par_all36_prot_lipid.prm   &nbsp;&nbsp;&nbsp;&nbsp;  #Force field file  <br/>
copy 200                                &nbsp;&nbsp;&nbsp;&nbsp;   #Copy number  <br/>

The method used here is to refer to the paper: <br/>
On Side-Chain Conformational Entropy of Proteins <br/>
