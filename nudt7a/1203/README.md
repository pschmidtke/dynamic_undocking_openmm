# How to set up
1. clone repo
2. put the module directory in your local PYTHONPATH environment variable
# How to run
1. Parametrize system:
  python 01_parametrize.py protein.pdb ligand.mol2
2. python 02_equil.py
3. python 03_md.py equil.chk md.chk md.csv md.pdb
4. mkdir smd_1
5. python 04_smd.py 300 equil.chk smd_1/out.csv smd_1/out.dat smd_1/out.pdb smd_1/out.dcd
