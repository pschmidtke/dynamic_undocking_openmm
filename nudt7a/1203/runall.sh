#!/bin/bash
source ~/.bash_profile

for i in {1..20}
  do
  	ni=$((i+1))
	python3 03_md.py md/md${i}.chk md/md${ni}.chk md/md${ni}.csv md/md${ni}.pdb 
	python3 04_smd.py 300 md/md${ni}.chk smd/out_300_{ni}.csv smd/out_300_{ni}.dat smd/out_300_{ni}.pdb smd/out_300_{ni}.dcd
	python3 04_smd.py 325 md/md${ni}.chk smd/out_300_{ni}.csv smd/out_300_{ni}.dat smd/out_300_{ni}.pdb smd/out_300_{ni}.dcd
  done