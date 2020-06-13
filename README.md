# VASP2ABINIT
a simple code to convert VASP POSCAR to ABINIT input file
The code was used in a project testing the time-reversal symmtery feature of LOBSTER as reported in this publication:
http://dx.doi.org/10.1002/jcc.26353
The program can also take an int argument to specify a kpoint density for automatic setup of the number of k-mesh.

The code needs the boost library and can be compiled by the following command:

g++ -o vasp2abinit.x vasp2abinit.cpp -lboost_filesystem  -lboost_system

Once the binary is generated, you need to put the binary in your bin directory. To run the program simply use the following command:

vasp2abinit.x \<path-to-poscar-or-contcar\>/POSCAR

OR 

vasp2abinit.x \<path-to-poscar-or-contcar\>/POSCAR DENSITY

where DENSITY is an integer number representing the kpoint density (in the unit of kpoints x atoms) to generate the number of kpoints for three directions.
The result will be output to the stdout. Once the ABINIT input file is produced, you still need to set the values for the variables ecut, pawecutdg, and nband, 
in order for the input file to be useable for a ABINIT calculation.

An example can be found in the directory "example". The output (mp-720.in) was produced with the following command:

vasp2abinit.x POSCAR_mp-720 8000 > mp-720.in
