import sys
import os
from ase import io, Atoms
from ase.data import atomic_numbers, atomic_masses
from ase.build import make_supercell, find_optimal_cell_shape, stack
from ase.constraints import FixAtoms, UnitCellFilter, ExpCellFilter, StrainFilter
from ase.optimize import QuasiNewton
from ase.calculators.emt import EMT
from numpy import *

basepath=os.getcwd()

pseudodir='/global/common/cori/software/vasp/pseudopotentials/PBE/potpaw_PBE.54/'

#infile1 = sys.argv[2]
#infile2 = sys.argv[3]

vdw_kernel_file='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/vdw_kernel.bindat'

xc_list = [ 'revvdWDF2', 'SCAN+rVV10', 'DFT-D3']

iop=sys.argv[1]

if iop == 'rotate' : rot_list = [ 1.0, 3.0, 5.0, 7.0, 10.0, 20.0 ]
#if iop == 'rotate_AB' : rot_list = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 20.0 ]
#if iop == 'rotate_AA' : rot_list = [ 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0, 190.0, 200.0 ]

def read_POSCAR(xc):
    infile=os.path.join(os.path.join(basepath,'bilayer'),'%s/CONTCAR')%(xc)
    infile_exists=os.path.exists(infile)
    if not infile_exists:
        print("POSCAR file does not exists")
    rf = open(infile,'rt')
    comment_line=rf.readline()
    lattice_scale=float(rf.readline().split()[0])
    lattice_vectors=[]
    for line in range(3):
        lat_abc=rf.readline().split()
        floatvec=float(lat_abc[0]),float(lat_abc[1]),float(lat_abc[2])
        lattice_vectors.append(floatvec)

    basis_vectors=array(lattice_vectors)*lattice_scale
    cellX=basis_vectors[0]
    cellY=basis_vectors[1]
    cellZ=basis_vectors[2]
    cell=(cellX, cellY, cellZ)

    atom_symbols=[]
    numofatoms=rf.readline().split()
    atomtypes=numofatoms
    numofatoms=rf.readline().split()
    numsys = len(numofatoms)
    for i, num in enumerate(numofatoms):
        numofatoms[i] = int(num)
        [atom_symbols.append(atomtypes[i]) for na in range(numofatoms[i])]
    coord_sys=rf.readline()
    cartesian = coord_sys[0].lower() == 'c' or coord_sys[0].upper() == 'C'
    tot_natoms=sum(numofatoms)
    atoms_pos=empty((tot_natoms,3))
    for atom in range(tot_natoms):
        ap = rf.readline().split()
        atoms_pos[atom]=(float(ap[0]),float(ap[1]),float(ap[2]))
    if cartesian:
        atoms_pos *= lattice_scale
    atoms = Atoms(symbols=atom_symbols, cell=basis_vectors, pbc=True)
    if cartesian:
        atoms.set_positions(atoms_pos)
    else:
        atoms.set_scaled_positions(atoms_pos)
    #print(atoms,lattice_scale,cellX,cellY,cellZ)
    return atoms,lattice_scale,cellX,cellY,cellZ

#writing k-points
def write_kpoints():
    ofkp=open("KPOINTS",'w')
    kpoints='''Automatic Mesh 
0 
G(M) 
12 12 6
0  0  0
'''
    ofkp.write(kpoints)

#writing INCAR files
def write_incar(xc):
    ofinc= open("INCAR",'w')
    incar_data='''!System settings
   PREC   =  accurate     low, normal, accurate
   ISTART =      0     job   : 0-new  1-cont  2-samecut
   ICHARG =      2     charge: 0-INIWAV 1-CHGCAR 2-atom 10-nonself
!   INIWAV =      1     0-jellium 1-rand
   ISYM   =      0

! Magnetic Moment 
     ISPIN  =      2     spin polarized? 1:no 2:yes
     MAGMOM = 2*0 4*0  ! local magnetic moment parallel to SAXIS
!    SAXIS =  0 0 1   ! quantization axis parallel to vector (x,y,z)

! Spin-orbit Coupling
!  LSORBIT = .TRUE.
!  ICHARG = 11      ! non selfconsistent run, read CHGCAR
!  LMAXMIX = 4      ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6
                           ! you need to set LMAXMIX already in the collinear calculation
!  NBANDS = 2 * number of bands of collinear run

!Electronic Relaxation
   ENCUT  =   450.0 eV
   !ENAUG  =  2700.0 eV
   NELM   =     200     max e SC
   !NELMIN =      10     min e SC
   !NELMDL =     -5     non-SC +:every -:first
   EDIFF  =   1.0e-5     stopping-criterion for ELM
   LREAL  =   .FALSE.   !Auto    real-space projection
   ALGO   =   fast     normal:38(D) fast:38-48(RMM) all:58 damped:53

   ISMEAR =     0     1,2:metal -5,0:sem/ins
   SIGMA  =   0.05     broadening in eV
'''
    if xc == "PBEsol" :
       incar_data +='''  
    GGA = PS 
'''
    elif xc == "optB88" :
        incar_data +=''' 
    GGA = BO
    PARAM1 = 0.1833333333
    PARAM2 = 0.2200000000
    LUSE_VDW = .TRUE.
    AGGAC = 0.0000
    LASPH = .TRUE.
'''
    elif xc == "revvdWDF2" :
        incar_data +='''
    GGA      = MK
    !GGA      = ML
    LUSE_VDW = .TRUE.
    PARAM1   = 0.1234
    PARAM2   = 0.711357
    Zab_vdW  = -1.8867
    AGGAC    = 0.0000
    LASPH = .TRUE.
'''
    elif xc == "SCAN+rVV10" :
        incar_data += '''
    METAGGA  = SCAN
    LUSE_VDW = .TRUE.
    BPARAM = 15.7     # default but can be overwritten by this tag
    CPARAM = 0.0093  # default but can be overwritten by this tag
    LASPH = .TRUE.
'''
    elif xc == "DFT-TS" :
        incar_data += '''
    IVDW= 2 | 20
    !GGA = PE
    !LVDWSCS=.TRUE.
    !ADDGRID=.FALSE.
   ! LSCSGRAD=.TRUE.
   ! LSCALER0=.TRUE.
'''
    elif xc == "DFT-D3" :
        incar_data +='''  
   IVDW=12
   GGA = PE
'''
    incar_data += '''
!Ionic relaxation
   EDIFFG =   1.0e-3     stopping-criterion for IOM +:energy -:force
   NSW    =     0        number of steps for IOM
   IBRION =    -1       ionic relax: 0-MD 1-static 2-CG 3-dampedMD
   !POTIM  =   0.5     step size for ionic-motion
   !ISIF   =   4        2:relax ions only; 3:also relax volume and cell shape; 4:relax ions+cellshape, volume=fixed
   !TEBEG  =  0.0
   !TEEND  =  0.0     temperature during run
!   SMASS  =     -1     -3:micro -2:const -1:scaled >=0:Nose mass


   LWAVE= .TRUE.
   LCHARG = .TRUE.
   NWRITE = 3

!For GPU; NPAR must be equal to #MPI_ranks(porcessors)/KPAR
!   LscaAWARE=.FALSE.
!   NPAR = 32
!   NSIM = 1

! Wannier 90 functionallity
! LWANNIER90=.TRUE. 
'''
    ofinc.write(incar_data)

#writing files needed for post-processing

#Preparing non self-consistent calculation
    
#a) with soc
    
def sc_soc():
    chng_sc_soc = open("change_sc_soc.sh",'w')
    chng_sc_soc_data='''#! /bin/bash/

cp -r nscf_wsoc/CONTCAR_nscf_wsoc POSCAR
cp -r nscf_wsoc/KPOINTS_nscf KPOINTS

python change_sc_input_soc.py

'''
    chng_sc_soc.write(chng_sc_soc_data)

#Changing input files for nscf calculations

#a) With SOC
def change_input_sc_soc():
    chng_sc_inpt_soc = open("change_sc_input_soc.py",'w')
    chng_sc_inpt_soc_data="""import sys
import os

in_file = open("nscf_wsoc/INCAR_nscf_wsoc", "rt")
out_file = open("INCAR", "wt")


for line in in_file:
    line=line.replace('MAGMOM = 2*0 4*0',  'MAGMOM = 0 0 2 0 0 -2 12*0')
    line=line.replace('!    SAXIS',  '     SAXIS')
    line=line.replace('!  LSORBIT',  '     LSORBIT')
    line=line.replace('ISIF','!ISIF')
    out_file.write(line)

in_file.close()
out_file.close()

out_file = open("INCAR", "a")

append_line='''
!DOS Calculations
    LORBIT = 11      output DOSCAR (DOS) and PROCAR (PDOS)
    LORBMOM =.TRUE.
!   NEDOS = 1000     number of points for DOS
!   EMIN = -5       boundaries for energy range for DOS
!   EMAX = 5
!   NBANDS= *

!SOC extra flags
   LMAXMIX = 4      ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6
                           ! you need to set LMAXMIX already in the collinear calculation
   GGA_COMPAT = .FALSE.
   ADDGRID = .TRUE.
!  NBANDS = 2 * number of bands of collinear run

'''
out_file.write(append_line)

out_file.close()
"""
    chng_sc_inpt_soc.write(chng_sc_inpt_soc_data)

#Preparing DOS calculation

#a) Without SOC
def dos_wsoc():
    chng_dos_wsoc=open("change_dos_wsoc.sh",'w')
    chng_dos_wsoc_data='''#! /bin/bash/

file=nscf_wsoc

[ ! -f $file ] && mkdir $file

cp -r INCAR nscf_wsoc/INCAR_nscf_wsoc
cp -r OUTCAR nscf_wsoc/OUTCAR_nscf_wsoc
#cp -r WAVECAR nscf_wsoc/WAVECAR_nscf_wsoc
#cp -r CHGCAR nscf_wsoc/CHGCAR_nscf_wsoc
#cp -r vasprun.xml nscf_wsoc/vasprun_nscf_wsoc.xml
cp -r KPOINTS nscf_wsoc/KPOINTS_nscf
cp -r CONTCAR nscf_wsoc/CONTCAR_nscf_wsoc
cp -r POSCAR nscf_wsoc/POSCAR_nscf_wsoc

#cp -r *.o nscf_wsoc/

python change_dos_input_wsoc.py

'''
    chng_dos_wsoc.write(chng_dos_wsoc_data)

#b) With SOC

def dos_soc():
    chng_dos_soc=open("change_dos_soc.sh",'w')
    chng_dos_soc_data='''#! /bin/bash/

file=nscf_soc

[ ! -f $file ] && mkdir $file

cp -r INCAR nscf_soc/INCAR_nscf_soc
cp -r OUTCAR nscf_soc/OUTCAR_nscf_soc
#cp -r WAVECAR nscf_soc/WAVECAR_nscf_soc
#cp -r CHGCAR nscf_soc/CHGCAR_nscf_soc
#cp -r vasprun.xml nscf_soc/vasprun_nscf_soc.xml
cp -r KPOINTS nscf_soc/KPOINTS_nscf
cp -r CONTCAR nscf_soc/CONTCAR_nscf_soc

#cp -r *.o nscf_soc/

python change_dos_input_soc.py

'''
    chng_dos_soc.write(chng_dos_soc_data)

#Preparing input files for DOS calculations

#a) Without SOC

def change_input_dos_wsoc():    
    chng_dos_input_wsoc=open("change_dos_input_wsoc.py",'w')
    chng_dos_input_wsoc_data='''import sys
import os

in_file = open("nscf_wsoc/INCAR_nscf_wsoc", "rt")
out_file = open("INCAR", "wt")


for line in in_file:
    line=line.replace('ISTART =      0','ISTART =      1')
    line=line.replace('ICHARG =      2','ICHARG =      11')
    line=line.replace('ISMEAR =     0 ','ISMEAR =     -5')
    line=line.replace('SIGMA  =   0.05 ','SIGMA  =   0.0')
    line=line.replace('LWAVE= .TRUE.','LWAVE= .FALSE.')
    out_file.write(line)

in_file.close()
out_file.close()

out_file = open("INCAR", "a")
'''
    chng_dos_input_wsoc_data +="""
append_line='''
!DOS Calculations
   LORBIT = 11      output DOSCAR (DOS) and PROCAR (PDOS)
   NEDOS = 1000     number of points for DOS
!   EMIN = -5       boundaries for energy range for DOS
!   EMAX = 5
!   NBANDS= *
'''
out_file.write(append_line)

out_file.close()

in_file_kpt = open("nscf_wsoc/KPOINTS_nscf", "rt")
out_file_kpt = open("KPOINTS", "wt")

for line in in_file_kpt:
    line=line.replace('12 12 6','25 25 6')
    out_file_kpt.write(line)

in_file_kpt.close()
out_file_kpt.close()
"""
    chng_dos_input_wsoc.write(chng_dos_input_wsoc_data)

#b) With SOC
    
def change_input_dos_soc():
    chng_dos_input_soc=open("change_dos_input_soc.py",'w')
    chng_dos_input_soc_data='''import sys
import os

in_file = open("nscf_soc/INCAR_nscf_soc", "rt")
out_file = open("INCAR", "wt")


for line in in_file:
    line=line.replace('ISTART =      0','ISTART =      1')
    line=line.replace('ICHARG =      2','ICHARG =      11')
    line=line.replace('ISMEAR =     0 ','ISMEAR =     -5')
    line=line.replace('SIGMA  =   0.05 ','SIGMA  =   0.0')
    line=line.replace('LWAVE= .TRUE.','LWAVE= .FALSE.')
    line=line.replace('!   NEDOS','    NEDOS')
    out_file.write(line)

in_file.close()
out_file.close()

in_file_kpt = open("nscf_soc/KPOINTS_nscf", "rt")
out_file_kpt = open("KPOINTS", "wt")

for line in in_file_kpt:
    line=line.replace('12 12 6','25 25 6')
    out_file_kpt.write(line)

in_file_kpt.close()
out_file_kpt.close()
'''
    chng_dos_input_soc.write(chng_dos_input_soc_data)

# Processing DOS 

#a) Without SOC

def process_dos_wsoc(atom1,atom2):
    proces_dos_wsoc=open("process_dos_wsoc.sh",'w')
    proces_dos_wsoc_data='''#! /bin/bash/

#echo -e "11\\n113\\n"|vaspkit
#echo -e "11\\n111\\n"|vaspkit

file=dos_wsoc

[ ! -f $file ] && mkdir $file

cp -r INCAR dos_wsoc/INCAR_dos_wsoc
cp -r OUTCAR dos_wsoc/OUTCAR_dos_wsoc
#cp -r vasprun.xml dos_wsoc/vasprun_dos_wsoc.xml
cp -r DOSCAR dos_wsoc/DOSCAR_wsoc
#cp -r PROCAR dos_wsoc/PROCAR_wsoc
cp -r KPOINTS dos_wsoc/KPOINTS_dos_wsoc
cp -r POSCAR  dos_wsoc/POSCAR_dos_wsoc

#mv *.o dos_wsoc/

#mv TDOS.dat dos_wsoc/TDOS_wsoc.dat
#mv ITDOS.dat dos_wsoc/ITDOS_wsoc.dat

#for i in {UP,DW}; do
#for j in {%s,%s}; do
#mv IPDOS_$j'_'$i.dat dos_wsoc/IPDOS_$j'_'$i'_wsoc'.dat
#mv PDOS_$j'_'$i.dat dos_wsoc/PDOS_$j'_'$i'_wsoc'.dat
#done
#done
'''%(atom1,atom2)

    proces_dos_wsoc.write(proces_dos_wsoc_data)

#b) With SOC

def process_dos_soc(atom1,atom2):    
    proces_dos_soc=open("process_dos_soc.sh",'w')
    proces_dos_soc_data='''#! /bin/bash/

#echo -e "11\\n113\\n"|vaspkit
#echo -e "11\\n111\\n"|vaspkit

file=dos_soc

[ ! -f $file ] && mkdir $file

cp -r INCAR dos_soc/INCAR_dos_soc
cp -r OUTCAR dos_soc/OUTCAR_dos_soc
#cp -r vasprun.xml dos_soc/vasprun_dos_soc.xml
cp -r DOSCAR dos_soc/DOSCAR_soc
#cp -r PROCAR dos_soc/PROCAR_soc
cp -r KPOINTS dos_soc/KPOINTS_dos_soc
cp -r POSCAR  dos_soc/POSCAR_dos_soc

#mv *.o dos_soc/

#mv TDOS.dat dos_soc/TDOS_soc.dat
#mv ITDOS.dat dos_soc/ITDOS_soc.dat
#mv I%s.dat dos_soc/I%s_soc.dat
#mv I%s.dat dos_soc/I%s_soc.dat
#mv *_SOC* dos_soc/
'''%(atom1,atom1,atom2,atom2)

    proces_dos_soc.write(proces_dos_soc_data)

# Preparing Band Structure Calculations

#a) Without SOC
   
def bs_wsoc():
    chng_bs_wsoc=open("change_band_str_wsoc.sh",'w')
    chng_bs_wsoc_data='''#! /bin/bash/

sed -i 's/ISMEAR =     -5/ISMEAR =     0/' INCAR
sed -i 's/SIGMA  =   0.0/SIGMA  =   0.05/' INCAR

echo -e "3\\n303\\n" |vaspkit

cp -r KPATH.in KPOINTS

sed -i 's/20/100/' KPOINTS
'''
    chng_bs_wsoc.write(chng_bs_wsoc_data)

#b) With SOC
    
def bs_soc():
    chng_bs_soc=open("change_band_str_soc.sh",'w')
    chng_bs_soc_data='''#! /bin/bash/

sed -i 's/ISMEAR =     -5/ISMEAR =     0/' INCAR
sed -i 's/SIGMA  =   0.0/SIGMA  =   0.05/' INCAR

echo -e "3\\n303\\n" |vaspkit

cp -r KPATH.in KPOINTS

sed -i 's/20/100/' KPOINTS
'''
    chng_bs_soc.write(chng_bs_soc_data)

# Processing Band Structure

#a) Without SOC

def process_bs_wsoc(atom1,atom2):
    proces_bs_wsoc=open("process_band_str_wsoc.sh",'w')
    proces_bs_wsoc_data = '''#! /bin/bash/

#echo -e "21\\n213\\n"|vaspkit
#echo -e "21\\n211\\n"|vaspkit

file=band_str_wsoc

[ ! -f $file ] && mkdir $file

cp -r INCAR band_str_wsoc/INCAR_bs_wsoc
cp -r OUTCAR band_str_wsoc/OUTCAR_bs_wsoc
cp -r KPOINTS band_str_wsoc/KPOINTS_bs_wsoc
cp -r POSCAR band_str_wsoc/POSCAR_bs_wsoc
cp -r DOSCAR band_str_wsoc/DOSCAR_bs_wsoc
cp -r PROCAR band_str_wsoc/PROCAR_bs_wsoc
cp -r EIGENVAL band_str_wsoc/EIGENVAL_bs_wsoc
cp -r vasprun.xml band_str_wsoc/vasprun_bs_wsoc.xml

#mv *.o band_str_wsoc/

#mv BAND.dat band_str_wsoc/BAND_wsoc.dat
#mv BAND_GAP band_str_wsoc/BAND_GAP_wsoc.dat
#mv KLABELS band_str_wsoc/KLABELS_wsoc
#mv KLINES.dat band_str_wsoc/KLINES_wsoc.dat

#for i in {UP,DW}; do
#mv REFORMATTED_BAND_$i.dat band_str_wsoc/REFORMATTED_BAND_$i'_wsoc'.dat
#for j in {%s,%s}; do
#mv PBAND_$j'_'$i.dat band_str_wsoc/PBAND_$j'_'$i'_wsoc'.dat
#done
#done
'''%(atom1,atom2)

    proces_bs_wsoc.write(proces_bs_wsoc_data)

#b) With SOC
    
def process_bs_soc():
    proces_bs_soc = open("process_band_str_soc.sh",'w')
    proces_bs_soc_data='''#! /bin/bash/

#echo -e "21\\n213\\n"|vaspkit
#echo -e "21\\n211\\n"|vaspkit

file=band_str_soc

[ ! -f $file ] && mkdir $file

cp -r INCAR band_str_soc/INCAR_bs_soc
cp -r OUTCAR band_str_soc/OUTCAR_bs_soc
cp -r KPOINTS band_str_soc/KPOINTS_bs_soc
cp -r POSCAR band_str_soc/POSCAR_bs_soc
cp -r DOSCAR band_str_soc/DOSCAR_bs_soc
cp -r PROCAR band_str_soc/PROCAR_bs_soc
cp -r EIGENVAL band_str_soc/EIGENVAL_bs_soc
cp -r vasprun.xml band_str_soc/vasprun_bs_soc.xml

#mv *.o band_str_soc/

#mv BAND.dat band_str_soc/BAND_soc.dat
#mv BAND_GAP band_str_soc/BAND_GAP_soc.dat
#mv KLABELS band_str_soc/KLABELS_soc
#mv KLINES.dat band_str_soc/KLINES_soc.dat
#mv REFORMATTED_BAND.dat band_str_soc/REFORMATTED_BAND_soc.dat
#mv *_SOC* band_str_soc/

cp -r nscf_wsoc/POSCAR_nscf_wsoc POSCAR
'''
    proces_bs_soc.write(proces_bs_soc_data)  

#Plotting DOS and Band Structure

def plot_dos_bs():
    plot_DOS_BS=open("plot_band_struc_dos.py",'w')
    plot_DOS_BS_data='''from pymatgen.io.vasp import Vasprun, BSVasprun
from pymatgen.electronic_structure.plotter import BSDOSPlotter
import sys,os

iop1 = sys.argv[1]
iop2 = sys.argv[2]
iop3 = sys.argv[3]

v = BSVasprun(iop1)
v1 = Vasprun(iop1)
cdos = v1.complete_dos
bs = v.get_band_structure(kpoints_filename=iop2,line_mode=True)
plot = BSDOSPlotter(bs_projection='elements', dos_projection='elements')
plt= plot.get_plot(bs,dos=cdos)
plt.savefig(iop3)
plt.show()
'''
    plot_DOS_BS.write(plot_DOS_BS_data)

#Preparing postprocessing for submit file

def submit_pp(atom1,atom2,iop,xc):
    submit_pp= open("submit_pp.sh",'w')
    submit_pp_data='''#!/bin/bash
#SBATCH --job-name="%s-%s-%s-%s"
#SBATCH --output="%s-%s-%s-%s.o"
#SBATCH --error="%s-%s-%s-%s.o"
#SBATCH --account="m2663"
#SBATCH --qos=regular
#SBATCH -N 4
#SBATCH --constraint=knl
#SBATCH --time=10:00:00
#SBATCH --mail-user="debajit.chakraborty@temple.edu"
#SBATCH --mail-type=FAIL
#
#SBATCH --comment=15:00:00
#SBATCH --time-min=10:0:00
#SBATCH --signal=B:USR1@300
#SBATCH --requeue
#SBATCH --open-mode=append

module load openmpi python vasp/20181030-knl
export OMP_NUM_THREADS=8

#Run non-scf calculation (without SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Prepare for the DOS (without SOC)
sh change_dos_wsoc.sh

#Run DOS calculation (without SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Process DOS files (without SOC)
sh process_dos_wsoc.sh

# Prepare for Band Structure Calcualation (without SOC)
sh change_band_str_wsoc.sh

#Run Band Structure calculation (without SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Process band structure files (without SOC)
sh process_band_str_wsoc.sh

#Plot the band-gap and DOS (without SOC)
python plot_band_struc_dos.py band_str_wsoc/vasprun_bs_wsoc.xml band_str_wsoc/KPOINTS_bs_wsoc %s-%s_%s_%s_wsoc.png 

#SOC calculation
# Prepare for the non self-consistent calculation (SOC)
sh change_sc_soc.sh

#Run non-scf calculation (SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_ncl

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Prepare for the DOS (SOC)
sh change_dos_soc.sh

#Run DOS calculation (SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_ncl

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Process DOS files (SOC)
sh process_dos_soc.sh

# Prepare for Band Structure Calculation (SOC)
sh change_band_str_soc.sh

#Run Band Structure calculation (SOC)
srun -n 32 -c 32  --cpu_bind=cores vasp_ncl

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

# Process band structure files (SOC)
sh process_band_str_soc.sh

#Plot the band-gap and DOS (SOC)
python plot_band_struc_dos.py band_str_soc/vasprun_bs_soc.xml band_str_soc/KPOINTS_bs_soc %s-%s_%s_%s_soc.png 
# requeueing the job if remaining time >0
. /usr/common/software/variable-time-job/setup.sh
requeue_job func_trap USR1

wait
''' %(atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc)

    submit_pp.write(submit_pp_data)

#Producing Wannier90 input_files

def wannier_90_wsoc():
    chng_wann_wsoc=open("change_wannier_wsoc.sh",'w')    
    chng_wann_wsoc_data='''#! /bin/bash/

cp -r nscf_wsoc/INCAR_nscf_wsoc INCAR
cp -r nscf_wsoc/KPOINTS_nscf KPOINTS
cp -r nscf_wsoc/CONTCAR_nscf_wsoc POSCAR

sed -i 's/ISPIN  =      2/ISPIN  =      0/' INCAR
sed -i 's/! LWANNIER90/ LWANNIER90/' INCAR 

echo -e "Begin Projections
	Random
End Projections

guiding_centres=true
" > wannier90.win

'''
    chng_wann_wsoc.write(chng_wann_wsoc_data)

def wannier_90_soc():
    chng_wann_soc=open("change_wannier_soc.sh",'w')    
    chng_wann_soc_data='''#! /bin/bash/

cp -r nscf_soc/INCAR_nscf_soc INCAR
cp -r nscf_soc/KPOINTS_nscf KPOINTS
cp -r nscf_soc/CONTCAR_nscf_soc POSCAR

sed -i 's/!LWANNIER90/LWANNIER90/' INCAR 

echo -e "Begin Projections
	Random
End Projections

guiding_centres=true
" > wannier90.win

'''
    chng_wann_soc.write(chng_wann_soc_data)

def save_wannier_90_wsoc(atom1,atom2,iop,xc):
    chng_wann_wsoc=open("save_wannier_wsoc.sh",'w')    
    chng_wann_wsoc_data='''#! /bin/bash/

file=wannier_wsoc

[ ! -f $file ] && mkdir $file

cp -r wannier90* $file/

cd $file/
mv wannier90.win %s%s-%s-%s.win
mv wannier90.eig %s%s-%s-%s.eig
mv wannier90.amn %s%s-%s-%s.amn
mv wannier90.mmn %s%s-%s-%s.mmn
mv wannier90.wout %s%s-%s-%s.wout
cd ../

'''
    chng_wann_wsoc.write(chng_wann_wsoc_data)

def save_wannier_90_soc(atom1,atom2,iop,xc):
    chng_wann_soc=open("save_wannier_soc.sh",'w')    
    chng_wann_soc_data='''#! /bin/bash/

file=wannier_soc

[ ! -f $file ] && mkdir $file

cp -r wannier* $file/

cd $file/
mv wannier90.win %s%s-%s-%s.win
mv wannier90.eig %s%s-%s-%s.eig
mv wannier90.amn %s%s-%s-%s.amn
mv wannier90.mmn %s%s-%s-%s.mmn
mv wannier90.wout %s%s-%s-%s.wout
cd ../

'''
    chng_wann_soc.write(chng_wann_soc_data)

def submit_wannier(atom1,atom2,iop,xc):
    submit_wann= open("submit_wann.sh",'w')
    submit_wann_data="""#!/bin/bash
#SBATCH --job-name="%s-%s-%s-%s"
#SBATCH --output="%s-%s-%s-%s.o"
#SBATCH --error="%s-%s-%s-%s.o"
#SBATCH --account="m2663"
#SBATCH --qos=regular
#SBATCH -N 4
#SBATCH --constraint=knl
#SBATCH --time=10:00:00
#SBATCH --mail-user="debajit.chakraborty@temple.edu"
#SBATCH --mail-type=FAIL
#
#SBATCH --comment=12:00:00
#SBATCH --time-min=10:0:00
#SBATCH --signal=B:USR1@300
#SBATCH --requeue
#SBATCH --open-mode=append

module load vasp/20181030-knl wannier
export OMP_NUM_THREADS=8

# Preparing for Wannier Without SOC
sh change_wannier_wsoc.sh

# Running the Wannier 

srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

#Saving files (without SOC)
sh save_wannier_wsoc.sh

cd wannier_wsoc/

run -n 32 -c 32 --cpu_bind=cores wannier90.x -pp %s%s-%s-%s.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

srun -n 32 -c 32 --cpu_bind=cores wannier90.x %s%s-%s-%s.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

cp -r %s%s-%s-%s.win %s%s-%s-%s_band.win

echo "
restart = plot

bands_plot = true

bands_plot_format = gnuplot

bands_num_points = 100

begin kpoint_path
   G    0.0000000000   0.0000000000   0.0000000000  C     0.3333333333   0.3333333333   0.0000000000           
  C_2  -0.3333333333   0.6666666667   0.0000000000  Y_2  -0.5000000000   0.5000000000   0.0000000000          
  Y_2  -0.5000000000   0.5000000000   0.0000000000  G     0.0000000000   0.0000000000   0.0000000000        
   G    0.0000000000   0.0000000000   0.0000000000  M_2  -0.5000000000   0.5000000000   0.5000000000
  M_2  -0.5000000000   0.5000000000   0.5000000000  D    -0.3333333333   0.6666666667   0.5000000000     
  D_2   0.3333333333   0.3333333333   0.5000000000  A     0.0000000000   0.0000000000   0.5000000000                  
   A    0.0000000000   0.0000000000   0.5000000000  G     0.0000000000   0.0000000000   0.0000000000        
  L_2   0.0000000000   0.5000000000   0.5000000000  G     0.0000000000   0.0000000000   0.0000000000          
   G    0.0000000000   0.0000000000   0.0000000000  V_2   0.0000000000   0.5000000000   0.0000000000           
end kpoint_path

"  >>  %s%s-%s-%s_band.win

cp -r  %s%s-%s-%s.chk %s%s-%s-%s_band.chk
cp -r  %s%s-%s-%s.eig %s%s-%s-%s_band.eig

srun -n 32 -c 32 --cpu_bind=cores wannier90.x %s%s-%s-%s_band.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

cd ../

# Preparing for Wannier SOC
sh change_wannier_soc.sh

# Running the Wannier

srun -n 32 -c 32  --cpu_bind=cores vasp_ncl

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

#Saving files (SOC)
sh save_wannier_soc.sh

cd wannier_soc/

run -n 32 -c 32 --cpu_bind=cores wannier90.x -pp %s%s-%s-%s.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

srun -n 32 -c 32 --cpu_bind=cores wannier90.x %s%s-%s-%s.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

cp -r %s%s-%s-%s.win %s%s-%s-%s_band.win

echo "
restart = plot

bands_plot = true

bands_plot_format = gnuplot

bands_num_points = 100

begin kpoint_path
   G    0.0000000000   0.0000000000   0.0000000000  C     0.3333333333   0.3333333333   0.0000000000           
  C_2  -0.3333333333   0.6666666667   0.0000000000  Y_2  -0.5000000000   0.5000000000   0.0000000000          
  Y_2  -0.5000000000   0.5000000000   0.0000000000  G     0.0000000000   0.0000000000   0.0000000000        
   G    0.0000000000   0.0000000000   0.0000000000  M_2  -0.5000000000   0.5000000000   0.5000000000
  M_2  -0.5000000000   0.5000000000   0.5000000000  D    -0.3333333333   0.6666666667   0.5000000000     
  D_2   0.3333333333   0.3333333333   0.5000000000  A     0.0000000000   0.0000000000   0.5000000000                  
   A    0.0000000000   0.0000000000   0.5000000000  G     0.0000000000   0.0000000000   0.0000000000        
  L_2   0.0000000000   0.5000000000   0.5000000000  G     0.0000000000   0.0000000000   0.0000000000          
   G    0.0000000000   0.0000000000   0.0000000000  V_2   0.0000000000   0.5000000000   0.0000000000           
end kpoint_path

"  >>  %s%s-%s-%s_band.win

cp -r  %s%s-%s-%s.chk %s%s-%s-%s_band.chk
cp -r  %s%s-%s-%s.eig %s%s-%s-%s_band.eig

srun -n 32 -c 32 --cpu_bind=cores wannier90.x %s%s-%s-%s_band.win

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

cd ../

# requeueing the job if remaining time >0
. /usr/common/software/variable-time-job/setup.sh
requeue_job func_trap USR1

wait

"""%(atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc)
    submit_wann.write(submit_wann_data)

def create_input_files(atoms,lattice_scale,cellX,cellY,cellZ,iop,xc,rotate):
    print(rotate)
    coord=atoms.get_positions()
    posX = coord[:,0]
    posY = coord[:,1]
    posZ = coord[:,2]
    old_cell=(cellX,cellY,cellZ)
    rotate*=pi/180
    rotation_matrix= array(((cos(rotate), sin(rotate), 0.),
                      (-sin(rotate), cos(rotate), 0.),
                      (0., 0., 1.)))
    new_cell= dot(rotation_matrix,old_cell)
    #print(old_cell, rotation_matrix, new_cell)

    #DeltaX = max(posX) - min(posX)
    #DeltaY = max(posY) - min(posY)
    
    #DeltaZ = max(posZ) - min(posZ)
    
    #cellx=sqrt(sum(cellX**2))
    #celly=sqrt(sum(cellY**2))
    #cellz=DeltaZ+20
    #cell= (cellX, cellY, cellZ)
    #cell= [cellx,celly,cellz,90, 90, 120]
    atoms.set_cell(new_cell)
    atoms.wrap()
    xmin=atoms.positions[:,0].min()
    xmax=atoms.positions[:,0].max()
    ymin=atoms.positions[:,1].min()
    ymax=atoms.positions[:,1].max()
    zmin=atoms.positions[:,2].min()
    zmax=atoms.positions[:,2].max()
    atoms.positions += (xmax-xmin, ymax-ymin, 0)
#    atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
#    ucf=UnitCellFilter(atoms)
#    strain=StrainFilter(atoms)
#    print (ucf, strain)
#    calc = EMT()
#    atoms.calc=calc
#    qn = QuasiNewton(ucf)
#    traj = Trajectory('MoS2.traj', 'w', atoms)
#    qn.attach(traj)
#    qn.run(fmax=0.05)

    ind_sort = argsort(atoms.get_chemical_symbols())
    symbols = array(atoms.get_chemical_symbols())[ind_sort]
    new_coord = coord[ind_sort]
    outfile="POSCAR"

#writing a description of the output file
    new_comment=("%s"+""+"unit rotated bilayer")%(rotate*180/pi)
    
    of = open(outfile,'w')
    of.write(new_comment + '\n')

#writing the basis vectors for the new cell

    long_format=True

    of.write('%19.16f\n' % lattice_scale)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms.get_cell():
        of.write(' ')
        for el in vec:
            of.write(latt_form % el)
        of.write('\n')

#writing the atom symbols
    sc = []
    psym = symbols[0]
    count = 0
    for sym in symbols:
        if sym != psym:
            sc.append ((psym,count))
            psym = sym
            count = 1
        else:
            count += 1
    sc.append((psym,count))

    for sym, _ in sc:
        of.write(' {:3s}'.format(sym))
    of.write('\n')

    for _, count in sc:
        of.write(' {:3d}'.format(count))
    of.write('\n')

#writing the type of cell used
    of.write('Cartesian\n')
#writing the new coordinates
    if long_format:
        cform = ' %19.16f'
    else:
        cform = ' %9.6f'
    count=0
    for iatom, atom in enumerate(new_coord):
        for dcoord in atom:
            of.write(cform % dcoord)
        of.write('\n')
        count += 1
    tot_natoms=count

#writing POTCAR files
    ofpot=open("POTCAR",'w')
    for sym, _ in sc:
        psfile=os.path.join(pseudodir,'%s/POTCAR')%(sym)
        ofpot.write(open(psfile).read())

    atom1=symbols[0]
    atom2=symbols[tot_natoms-1]

    write_kpoints() 

    write_incar(xc)

    sc_soc()
    
    change_input_sc_soc()
   
    change_input_dos_wsoc()

    dos_wsoc()

    change_input_dos_soc()

    dos_soc()

    process_dos_wsoc(atom1,atom2)

    process_dos_soc(atom1,atom2)

    bs_wsoc()

    bs_soc()

    process_bs_wsoc(atom1,atom2)

    process_bs_soc()

    plot_dos_bs()

    submit_pp(atom1,atom2,iop,xc)

    wannier_90_wsoc()
    
    save_wannier_90_wsoc(atom1,atom2,iop,xc)
    
    wannier_90_soc()
    
    save_wannier_90_soc(atom1,atom2,iop,xc)

    submit_wannier(atom1,atom2,iop,xc)



new_dir = os.path.join(basepath,'%s')%(iop)
exist = os.path.exists(new_dir)
if not exist: os.mkdir(new_dir)
os.chdir(new_dir)
for k in xc_list :
    new_dir_1= os.path.join(new_dir,'%s')%(k)
    exist_1=os.path.exists(new_dir_1)
    if not exist_1:
        os.mkdir(new_dir_1)
    os.chdir(new_dir_1)
    atoms1,lat_scale,cell1X,cell1Y,cell1Z=read_POSCAR(k)
    print(atoms1.get_positions())
    for i in rot_list:
        rotate= i
        new_dir_2= os.path.join(new_dir_1,'%s_unit_rotated')%(i)
        exist_2=os.path.exists(new_dir_2)
        if not exist_2:
            os.mkdir(new_dir_2)
        os.chdir(new_dir_2)
        atoms11=atoms1[0,2,3]
        atoms12=atoms1[1,4,5]
        atoms11.get_positions()[:,0]=atoms11.get_positions()[:,0]-0
        atoms11.get_positions()[:,1]=atoms11.get_positions()[:,1]-0
        atoms12.get_positions()[:,0]=atoms12.get_positions()[:,0]-0
        atoms12.get_positions()[:,1]=atoms12.get_positions()[:,1]-0
        for j in range(len(atoms12)):
            print(atoms12.get_positions()[j])
        point=atoms12.get_positions()[0]
        atoms12.rotate(rotate,'z')
        #atoms12.euler_rotate(phi=rotate,center=point)
        print(atoms12.get_positions())
        atoms_rotated=atoms11+atoms12
        #atoms_rotated.wrap()
        create_input_files(atoms_rotated,lat_scale,cell1X,cell1Y,cell1Z,iop,k,i)
        if k == 'optB88' or k == 'revvdWDF2' : 
            vdW_file_dst= os.path.join(new_dir_2,'vdw_kernel.bindat')
            exist_3=os.path.exists(vdW_file_dst)
            if not exist_3:
                os.symlink(vdw_kernel_file,vdW_file_dst)
        print(os.getcwd())
        print (iop,k,i)
    print(os.getcwd())
print(os.getcwd())
