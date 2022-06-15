import sys
import os
from ase import io, Atoms
from ase.data import atomic_numbers, atomic_masses
from numpy import *

basepath=os.getcwd()

pseudodir='/global/common/cori/software/vasp/pseudopotentials/PBE/potpaw_PBE.54/'

infile1 = sys.argv[2]
#infile2 = sys.argv[3]

vdw_kernel_file='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/vdw_kernel.bindat'

#xc_list = [ 'PBEsol', 'optB88', 'revvdWDF2', 'SCAN+rVV10', 'DFT-TS', 'DFT-D3']
xc_list = [ 'revvdWDF2', 'SCAN+rVV10', 'DFT-D3']

iop=sys.argv[1]

def read_POSCAR(infile):
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
12 12 12 
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
!   ISYM   =      0

! Magnetic Moment 
     ISPIN  =      2     spin polarized? 1:no 2:yes
     MAGMOM = 3*0 6*0  ! local magnetic moment parallel to SAXIS
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
   NSW    =    200        number of steps for IOM
   IBRION =    2       ionic relax: 0-MD 1-static 2-CG 3-dampedMD
   !POTIM  =   0.5     step size for ionic-motion
   ISIF   =   4        2:relax ions only; 3:also relax volume and cell shape; 4:relax ions+cellshape, volume=fixed
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

#Relaxing the volume
def change_incar_isif3(): 
    chng_isif3 = open("change_incar_isif3.sh",'w')
    chng_isif3_data='''#! /bin/bash/

file=relax_fv

[ ! -f $file ] && mkdir $file

cp -r INCAR relax_fv/INCAR_relax_fv

mv OUTCAR relax_fv/OUTCAR_relax_fv

#mv WAVECAR relax_fv/WAVECAR_relax_fv
#mv CHGCAR relax_fv/CHGCAR_relax_fv
#mv vasprun.xml relax_fv/vasprun_relax_fv.xml

#cp -r *.o relax_fv/

cp -r POSCAR relax_fv/POSCAR_init
cp -r CONTCAR relax_fv/CONTCAR_relax_fv
cp -r CONTCAR POSCAR

sed -i 's/ISIF   =   4/ISIF   =   3/' INCAR

'''
    chng_isif3.write(chng_isif3_data)


#Relaxing the ions only

#Preparing input_files

def change_incar_ion_relax(): 
    chng_ion_relax = open("change_incar_irelax.sh",'w')
    chng_ion_relax_data='''#! /bin/bash/

file=relax_vol

[ ! -f $file ] && mkdir $file

cp -r INCAR relax_vol/INCAR_relax_vol

mv OUTCAR relax_vol/OUTCAR_relax_vol

#mv WAVECAR relax_vol/WAVECAR_relax_vol
#mv CHGCAR relax_vol/CHGCAR_relax_vol
#mv vasprun.xml relax_fv/vasprun_relax_vol.xml

#cp -r *.o relax_vol/

cp -r POSCAR relax_vol/POSCAR
cp -r CONTCAR relax_vol/CONTCAR_relax_vol
cp -r CONTCAR POSCAR

sed -i 's/ISIF   =   3/ISIF   =   2/' INCAR

'''
    chng_ion_relax.write(chng_ion_relax_data)

#writing the submit files

def submit_optimze(atom1,atom2,iop,xc):
    submit= open("submit.sh",'w')
    submit_data='''#!/bin/bash
#SBATCH --job-name="%s-%s-%s-%s"
#SBATCH --output="%s-%s-%s-%s.o"
#SBATCH --error="%s-%s-%s-%s.o"
#SBATCH --account="m2663"
#SBATCH --qos=regular
#SBATCH -C knl
#SBATCH -N 4
#SBATCH --time=10:00:00
#SBATCH --mail-user="debajit.chakraborty@temple.edu"
#SBATCH --mail-type=FAIL
#
#SBATCH --comment=12:00:00
#SBATCH --time-min=10:0:00
#SBATCH --signal=B:USR1@300
#SBATCH --requeue
#SBATCH --open-mode=append

module load vasp/20181030-knl
export OMP_NUM_THREADS=8

# Running the fixed volume relaxation

srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

#Running the volume relaxation

sh change_incar_isif3.sh

srun -n 32 -c 32  --cpu_bind=cores vasp_std

#wait until VASP to complete
srun_pid=`ps -fle|grep srun|head -1|awk '{print $4}'`
echo srun pid is $srun_pid  >&2
wait $srun_pid

#Running ion relaxation
sh change_incar_irelax.sh

srun -n 32 -c 32  --cpu_bind=cores vasp_std

# requeueing the job if remaining time >0
. /usr/common/software/variable-time-job/setup.sh
requeue_job func_trap USR1

wait

''' %(atom1,atom2,iop,xc,atom1,atom2,iop,xc,atom1,atom2,iop,xc)
    submit.write(submit_data)


def create_input_files(atoms,lattice_scale,cellX,cellY,cellZ,iop,xc):
    coord=atoms.get_positions()
    #posX = coord[:,0]
    #posY = coord[:,1]
    #posZ = coord[:,2]
    #DeltaX = max(posX) - min(posX)
    #DeltaY = max(posY) - min(posY)
    #DeltaZ = max(posZ) - min(posZ)
    #cellZ[2]=DeltaZ+20
    cell= (cellX, cellY, cellZ)
    atoms.set_cell(cell)
    
    ind_sort = argsort(atoms.get_chemical_symbols())
    symbols = array(atoms.get_chemical_symbols())[ind_sort]
    new_coord = coord[ind_sort]
    outfile="POSCAR"

#writing a description of the output file
    new_comment="optimize structure"
    
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

    write_kpoints() 

    write_incar(xc)

    change_incar_isif3()

    change_incar_ion_relax()

    submit_optimze(symbols[0],symbols[tot_natoms-1],iop,xc)


atoms1,lat_scale,cell1X,cell1Y,cell1Z=read_POSCAR(infile1)

atoms11=atoms1[0:5:2]
atoms22=atoms1[1:6:2]
Delta=atoms22.get_center_of_mass()-atoms11.get_center_of_mass()
h=sqrt(sum(Delta**2))
atoms33=atoms11.copy()
atoms33.center(about=atoms11.positions[0])
count=0
for i in atoms33:
    atoms33.positions[count][2] += 2*h
    count +=1
print(atoms33.get_positions(), Delta, h)
atoms_shifted=atoms11+atoms22+atoms33
cell1Z[2] +=1.5*h

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
    create_input_files(atoms_shifted,lat_scale,cell1X,cell1Y,cell1Z,iop,k)
    if k == 'optB88' or k == 'revvdWDF2' : 
        vdW_file_dst= os.path.join(new_dir_1,'vdw_kernel.bindat')
        exist_3=os.path.exists(vdW_file_dst)
        if not exist_3:
            os.symlink(vdw_kernel_file,vdW_file_dst)
    print (iop,k)
    print(os.getcwd())
print(os.getcwd())

