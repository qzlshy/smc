# smc
Calculate protein side chain conformational entropy for a given protein backbone conformation.<br/>

Compile<br/>
make clean<br/>
make<br/>

Run<br/>
./long cal.conf<br/>
<br/>
cal.conf is the configuration file, defind the input files of the program. <br/>

 表头 | 表头  
-- | --
rotamer ./usefile/rotamer.pdb            | #rotamer library  
pdb ./ubq/ubq-1.pdb                      | #input protein structure, you can change to yourself file
res ./usefile/res.txt                    | #defind side chain dihedral angle of 20 type residues  
top ./usefile/top_all36_prot_lipid.rtf   | #Topology file  
prm ./usefile/par_all36_prot_lipid.prm   | #Force field file  
copy 200                                 | #Copy number  

The method used here is to refer to the paper: <br/>
On Side-Chain Conformational Entropy of Proteins <br/>
