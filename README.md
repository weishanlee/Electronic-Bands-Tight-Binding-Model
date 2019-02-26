# Electronic-Bands-Tight-Binding-Model
The codes in the repository intend to reproduce the classic work of electronic band structures of zincblende and diamond crystals by P. Vogl et al[1].

Descriptions of the fortran codes:

(1)The main program is tbhf_main.for.
(2)The tbhf.for calculates the eigenvalues and eigenvectors at K-points of the LCAO matrix whose elements are given in the tbparameters.dat file for different  materials
(3)The order of parameters, extracted from TABLE 1 in Ref [1], in the matrix elements is as follows:
    	E(s,a), E(p,a), E(s,c), E(p,c), V(s,s), V(x,x)
    	V(x,y), V(sa,pc), V(sc,pa), E(s*,a), E(s*,c)
    	V(s*a,pc), V(pa,s*c)
(4)Some regular eispack packages are required, including
	ch.for; 
	epslon.for; 
	htribk.for; 
	htridi.for; 
	pythang.for; 
	tql2.for; and
	tqlrat.for, which may be found via the website: http://www.netlib.org/eispack/. For copyright reasons the above fortran files are NOT uploaded and are NOT included in this repository.
(5) The selected k points are given in nkl_801.dat file. The electronic band strucure is given in eb.dat, while the density of state is in DOS.dat. 

The executable file is also included. To calculate other materials, one may try to change the 13 parameters in the tbparameters.dat file by just coping other parameters in the lists and pasting to replace the original parameters. Later, run the exe file. 

References:

[1]P. Vogl, H. P. Hjalmarson, and John D. Dow, "A Semi-empirical Tight-Binging Theory of the Electronic Structure of Semiconductors", J. Phys. Chem. Solids, Vol. 44, No.5 (1983) pp. 365-378.

