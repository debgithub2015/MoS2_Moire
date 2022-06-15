import sys
import os
from ase import io, Atoms
from ase.data import atomic_numbers, atomic_masses
from numpy import *

basepath='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/'

pseudodir='/global/common/cori/software/vasp/pseudopotentials/PBE/potpaw_PBE/'

xc_list = ['DFT-D3']

#xc_list = [ 'PBEsol', 'optB88', 'revvdWDF2', 'SCAN+rVV10', 'DFT-TS', 'DFT-D2', 'DFT-D3']

iop=sys.argv[1]

if iop == 'shift':
    sep_list = [ 0.5, 1.0, 1.5, 2.0, 3.0, 5.0] 
elif iop == 'rotate':
    rot_list = [ 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0 ]

def create_input_files(j,iop,xc):
    infile='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/POSCAR_MoS2_cart'
    #infile='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/bilayer/DFT-D3/bl/POSCAR'
    rf=open(infile)
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

    print(atoms.get_positions())

    atoms1 = atoms[0:5:2]
    atoms2 = atoms[1:6:2]
    

    if iop == 'shift':
        shift= j
        Delta = atoms2.get_center_of_mass() - atoms1.get_center_of_mass()
        dist = sqrt(sum(Delta**2))
        unit = Delta/dist
        print(atoms2.get_center_of_mass(),atoms1.get_center_of_mass(), Delta, dist, unit, shift*unit)
        atoms2.translate(shift*unit)
        print(atoms2.get_positions())
        atoms_shifted=atoms1+atoms2
    elif iop == 'rotate':
        rotation = j
        atoms2.euler_rotate(phi=rotation, center='COP')
        atoms_shifted=atoms1+atoms2
        print(atoms_shifted.get_positions())
        
    new_coord=atoms_shifted.get_positions()

#    posX = new_coord[:,0]
#    posY = new_coord[:,1]
#    posZ = new_coord[:,2]
#    #DeltaX = max(posX) - min(posX)
#    #DeltaY = max(posY) - min(posY)
#    DeltaZ = max(posZ) - min(posZ)
#    
#    cellZ[2]=DeltaZ+20
    cell= (cellX, cellY, cellZ)
    atoms_shifted.set_cell(cell)
    
    ind_sort = argsort(atoms_shifted.get_chemical_symbols())
    symbols = array(atoms_shifted.get_chemical_symbols())[ind_sort]
    new_coord = new_coord[ind_sort]
    #outfile="POSCAR"+"_"+str(j)+"_shifted"
    outfile="POSCAR"

#writing a description of the output file
    if iop == 'shift' :
        new_comment="structure is"+ '\t'+ str(j) +'\t'+ "unit shifted"
    elif iop == 'rotate' :
        new_comment="structure is"+ '\t'+ str(j) +'\t'+ "degree rotated"

    of = open(outfile,'w')
    of.write(new_comment + '\n')

#writing the basis vectors for the new cell

    long_format=True

    of.write('%19.16f\n' % lattice_scale)
    if long_format:
        latt_form = ' %21.16f'
    else:
        latt_form = ' %11.6f'
    for vec in atoms_shifted.get_cell():
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
    for iatom, atom in enumerate(new_coord):
        for dcoord in atom:
            of.write(cform % dcoord)
        of.write('\n')

#writing POTCAR files
    ofpot=open("POTCAR",'w')
    for sym, _ in sc:
        psfile=os.path.join(pseudodir,'%s/POTCAR')%(sym)
        ofpot.write(open(psfile).read())

#writing KPOINTS files
    ofkp=open("KPOINTS",'w')
    kpoints='''Automatic Mesh 
0 
G(M) 
12 12 6 
0  0  0
'''

    ofkp.write(kpoints)
 
#writing INCAR files
    ofinc= open("INCAR",'w')
    incar_data='''!System settings
   PREC   =  accurate     low, normal, accurate
   ISTART =      0     job   : 0-new  1-cont  2-samecut
   ICHARG =      2     charge: 0-INIWAV 1-CHGCAR 2-atom 10-nonself
   ISPIN  =      2     spin polarized? 1:no 2:yes
!   INIWAV =      1     0-jellium 1-rand
!   ISYM   =      0
! Magnetic Moment 
    MAGMOM = 2*0 4*0  ! local magnetic moment parallel to SAXIS
!    SAXIS = 0 0 1   ! quantization axis parallel to vector (x,y,z)

! Spin-orbit Coupling
!  LSORBIT = .TRUE.
!  ICHARG = 11      ! non selfconsistent run, read CHGCAR
!  LMAXMIX = 4      ! for d-elements increase LMAXMIX to 4, f-elements: LMAXMIX = 6
                           ! you need to set LMAXMIX already in the collinear calculation
!  SAXIS =  x y z   ! direction of the magnetic field
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
    BPARAM = 6.3     # default but can be overwritten by this tag
    CPARAM = 0.0093  # default but can be overwritten by this tag
    LASPH = .TRUE.
'''
    elif xc == "DFT-TS" :
        incar_data += '''
    IVDW=2
    GGA = PE
   ! LVDWSCS=.TRUE.
   ! LSCSGRAD=.TRUE.
   ! LSCALER0=.TRUE.
'''
    elif xc == "DFT-D2" :
        incar_data +='''  
   IVDW=1
   GGA = PE
'''
    elif xc == "DFT-D3" :
        incar_data +='''  
   IVDW=12
   GGA = PE
'''
    incar_data += '''
!Ionic relaxation
   EDIFFG =   1.0e-5     stopping-criterion for IOM +:energy -:force
   NSW    =     0        number of steps for IOM
   IBRION =    -1       ionic relax: 0-MD 1-static 2-CG 3-dampedMD
   !POTIM  =   0.5     step size for ionic-motion
   !ISIF   =   0     relax atomic coords only
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
'''
    ofinc.write(incar_data)

#writing the submit files
    submit= open("submit.sh",'w')
    submit_data='''#!/bin/bash
#SBATCH --job-name="%s"
#SBATCH --output="%s_%i.o"
#SBATCH --error="%s_%i.o"
#SBATCH --account="m2663"
#SBATCH --qos=regular
#SBATCH -N 4
#SBATCH --constraint=knl
#SBATCH --time=00:30:00
#SBATCH --mail-user="debajit.chakraborty@temple.edu"
#SBATCH --mail-type=FAIL

module load vasp/20181030-knl
export OMP_NUM_THREADS=8

srun -n 32 -c 32  --cpu_bind=cores vasp_std

'''%(iop,iop,j,iop,j)
    submit.write(submit_data)

if iop == 'shift' :
    new_dir = os.path.join(basepath,'%s')%(iop)
    exist = os.path.exists(new_dir)
    if not exist: os.mkdir(new_dir)
    os.chdir(new_dir)
    for k in xc_list :
        new_dir_1= os.path.join(new_dir,'%s')%(k)
        exist_1=os.path.exists(new_dir_1)
        if not exist_1: os.mkdir(new_dir_1)
        os.chdir(new_dir_1)
        for i in sep_list:
            new_dir_2=os.path.join(new_dir_1,'%s'+'_unit_shifted')%(i)
            exist_2= os.path.exists(new_dir_2)
            if not exist_2: os.mkdir(new_dir_2)
            os.chdir(new_dir_2)
            create_input_files(i,iop,k)
            print (i,iop,k)
            print(os.getcwd())
        print(os.getcwd())
    print(os.getcwd())
elif iop == 'rotate' :
    new_dir = os.path.join(basepath,'%s')%(iop)
    exist = os.path.exists(new_dir)
    if not exist: os.mkdir(new_dir)
    os.chdir(new_dir)
    for k in xc_list :
        new_dir_1= os.path.join(new_dir,'%s')%(k)
        exist_1=os.path.exists(new_dir_1)
        if not exist_1: os.mkdir(new_dir_1)
        os.chdir(new_dir_1)
        for i in rot_list:
            new_dir_2=os.path.join(new_dir_1,'%s'+'_deg_rot')%(i)
            exist_2= os.path.exists(new_dir_2)
            if not exist_2: os.mkdir(new_dir_2)
            os.chdir(new_dir_2)
            create_input_files(i,iop,k)
            print (i,iop,k)
            print(os.getcwd())
        print(os.getcwd())
    print(os.getcwd())
