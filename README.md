How to use the codes in mock_data_generation/v2/
================================================

Copy the folder mock_data_generation/v"$version"/template into a new directory e.g./ working, e.g.,

cd mock_data_generation/v"$version"/
cp -r template working
cd working

Proceed according to the instruction below.



First Part: Preparation of 'experimental.template.dat'
======================================================
A necessary template file to be used in the Second Part of mock data generation (see below) has to be first prepared based on a template file 'realistic_experimental.dat'. 'realistic_experimental.dat' is an output file obtained from an realistic FTIR experimental measurement. It contains a very lengthy list of wavelength from 400 to 1000 nm at about 1.3 nm interval. The total number of wavelength is 570. 

For the purpose of preparing mock data for GMT, the number of wavelength has to be trimmed to a smaller size, e.g. 57. To this end, read the instruction in trim_realistic_experimental.py file to edit the variable nskip. It is set to 10 by default to produce 57 wavelengths. Execute 

	python trim_realistic_experimental.py

to obtain the output file, e.g., experimental_trimmed.57.dat.

The file experimental_trimmed.57.dat will then be automatically copied and saved as 'experimental.template.dat'. 

Note that this Part provides a flexible way to set the resolution of the wavelengths to be used in the mock data generation in the Second Part, namely, experiental.dat. The resolution in the wavelengths in experiental.dat in the Second Part is an exact copy of the wavelength as recorded in 'experimental.template.dat'.

Proceed to the Second Part after you have prepared the 'experimental.template.dat' file. 


Second Part:
============

Mandatory requirement: Prexisting presence of input data file 'experimental.template.dat'. 
Three columns in 'experimental.template.dat' will be modified, namely,

column[1]		column[2]			column[3]			
Ref Index (real)	Ref Index (imaginary)		Cext

The Ref Index (real) and Ref Index (imaginary) columns are the real and imaginary parts of the reflective index of the individual sphere. In the simplest case, all spheres are of the same radius and same reflective index. In the original work by Xu Yi Lin, the originator of the GMM code, these spheres are made up of acryllic spheres known as 'bk7' with size expressed in units that is the same as that for the incident wavelength. The incident wavelengths that will be used are covered for a range as stated in the 'wavelength' column (i.e., column[0]) in experimental.dat. Note that the purpose of the file experimental.dat is merely to be used as a template for fixing the range and value of wavelength to be used when calculating the scattering coefficients by the GMT software. 

At the end of the following processes, a resultant file will be produced, 'experimental.dat', which has exactly the same data format as 'experimental.template.dat' but with the above-mentioned columns replaced by new data. In 'experimental.dat', the wavelength column (column[0]) is exactly the same as that in 'experimental.template.dat'. 

Step 1
Read hcp.py. Edit hcp.py for 
	noa     ## default 14
	lc      ## default 200 nm, commond diameter of the spheres
	n       ## detaul 1.5   Real part of reflective index
	nk      ## default 10.5 Imaginary  part of reflective index


Step 2
Execute 

	python hcp.py 
	
to generate random hcp configuration in the form of config.xyz and config.dat.
The spheres produced in config.dat has a common radius size, rc=lc/2. 
Input required in hcp.py: (1) number of atoms, noa and (2) lattice constatn, lc. 

Step 3
Execute 

	python create_experimental.dat.py 
	
to create experimental.dat from experimental.template.dat. 

At this stage, in the file 'experimental.dat', the Ref Index (real), Ref Index (imaginary) columns are obtained by replacing the original columns in 'experimental.template.dat' by the values as recorded in config.dat (which in turn was generated in hcp.py in Step 1). The file experimental.dat so created contains only the columns of wavelength, real refractive index and imaginary refractice index. The values of Cext column still remain as a string 'Cext'. They will be replaced by GMT-generated values only in Step 6.

Step 4
With the presence of config.dat produced in Step 2, execute convert_2_mycal.in.sh via

	python convert_config.dat2mycal.in.py

to produce radius.dat, mycal.in
	
Step 5
In the presence of mycal.in, experimental.dat, execute 

	./mycal.sh

As a result, gmt.dat will be produced. The gmt.dat so produced is the Cext vs. wavelength curve. The wavelengths are that listed in the first column in experimental.dat. 

Step 6
In the presence of gmt.dat obtained from the previous step, execute 

	python gen_mock_exp.py 

to produce an overwritten version of experimental.dat. At this stage, the Cext column (column[3]) from Step 3 is replaced by that from gmt.dat obtained from step 5.

Step 7
In the presence of plt.gmt_mock.norm.py, execute 

	python plt.gmt_mock.norm.py

to produce gmt_mock.norm.pdf	


Step 8. The process from Step 1 - Step 7 can be conveniently executed via an all-in-one script,

	./autorun.sh

The detailed input and output of the scripts described above are also described in autorun.sh. Execution of autorun.sh will produce a folder named in the form of e.g., 

	nwl.57.noa.14.lc.200.n.1.5.nk.10.5.0823160151

where the numbers are the specification of the variables nwl (number of wavelength), noa (number of spheres), lc (commond diameter of the spheres), n(real part of reflective index of the spheres), nk(imaginary part of the reflective index of the spheres). The figures 0823160151 is to be interpreted as follows:

	0823160151 = 08 (month=August) + 23(day) + 160151 (time in the format 16:01:51)

The folder such as nwl.57.noa.14.lc.200.n.1.5.nk.10.5.0823160151 will replace the name of the folder in which the autorun.sh package is executed. 
	

As a summary: Edit hcp.py to specify the number of spheres, their common radius, the real and imaginary parts of their common reflective index. With the presence of experimental.template.dat (which is prepared in Part B to fix the number of wavelengths), execute ./autorun.sh. The following output files will be produced in an new folder with indicative name such as nwl.57.noa.14.lc.200.n.1.5.nk.10.5.0823160151, namely:

	config.dat 
	config.xyz 
	radius.dat
	gmt_mock.dat
	gmt_mock.norm.dat
	experimental_mock.dat 
	experimental.dat
	gmt_mock.norm.pdf

These are the mock data required for dakota-gmt. Place these files into the directory common_fixed_radius/templatedir/. Continue to run dakota-gmt from there.
