#! /bin/bash 

#./clean.sh
rm -rf nwl.*
python hcp.py                          	## produce config.dat, config.xyz, specification statement

python create_experimental.dat.py      	## read in experimental.template.dat to produce experimental.dat
python convert_config.dat2mycal.in.py  	## read in config.dat, produce radius.dat, mycal.in
./mycal.sh    						   	## read in experimental.dat, produce gmt_mock.dat
python gen_mock_exp.py  				## read in gmt_mock.dat, experimental.dat (output from  	
										## create_experimental.dat.py), to produce 		
										## gmt_mock.norm.dat, experimental_mock.dat and experimental.dat. ## Here, the output experimental.dat (with Cext updated) is different than that produced by create_experimental.dat.py (in which the Cext column remains as 'Cext' string)
										## Also, experimental.dat=experimental_mock.dat

rm -rf gmt_mock.norm.pdf
python plt.gmt_mock.norm.py   			## read in gmt_mock.norm.dat, radius.dat, config.dat, produce 
										## gmt_mock.norm.pdf

dirname=$(ls nwl.* | awk 'END {print}')
echo $dirname
pwd=$(pwd)
echo $pwd

cwd=$pwd
cd ../
upath=$(pwd)
echo 'upath' $upath

name="$upath"/"$dirname"
echo $name
mv $pwd $name
#mv $pwd $dirname
