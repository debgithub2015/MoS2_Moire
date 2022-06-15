from ase import io, Atoms
from numpy import *
import sys

infile='/global/cscratch1/sd/debajitc/calculations/tmds/mos2/POSCAR_MoS2_cart'


#Reading the vasp POSCAR file
rf=open(infile)
comment_line=rf.readline()
lattice_scale=float(rf.readline().split()[0])
lattice_vectors=[]
for line in range(3):
    lat_abc=rf.readline().split()
    floatvec=float(lat_abc[0]),float(lat_abc[1]),float(lat_abc[2])
    lattice_vectors.append(floatvec)

basis_vectors=array(lattice_vectors)*lattice_scale

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
print(cartesian)
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
print(atoms)

#shifting units
sep_list = [ 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0]

cellX=basis_vectors[0]
cellY=basis_vectors[1]
cellZ=basis_vectors[2]


def create_POSCAR(j):
    shift= sep_list[j]
    
    atoms1 = atoms[0:5:2]
    atoms2 = atoms[1:6:2]

    print(atoms1,atoms2)

    Delta = atoms2.get_center_of_mass() - atoms1.get_center_of_mass()
    dist = sqrt(sum(Delta**2))
    unit = Delta/dist
    atoms2.translate(shift*unit)

    atoms_shifted=atoms1+atoms2
    new_coord=atoms_shifted.get_positions()
    print (atoms_shifted,new_coord)
    posX = new_coord[:,0]
    posY = new_coord[:,1]
    posZ = new_coord[:,2]
    #DeltaX = max(posX) - min(posX)
    #DeltaY = max(posY) - min(posY)
    DeltaZ = max(posZ) - min(posZ)
    cellZ[2]=DeltaZ+20
    cell= (cellX, cellY, cellZ)
    atoms_shifted.set_cell(cell)
    
    print(atoms_shifted)

    sort=True

    if sort:
        ind = argsort(atoms_shifted.get_chemical_symbols())
        symbols = array(atoms_shifted.get_chemical_symbols())[ind]
        new_coord = new_coord[ind]
    print(ind,symbols,new_coord)

    outfile="POSCAR"+"_"+str(j)+"_shifted"

#writing a description of the output file
    new_comment="structure is"+ '\t'+ str(j) +'\t'+ "unit shifted"

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
    print(sc)

    for sym, _ in sc:
        of.write(' {:3s}'.format(sym))
    of.write('\n')

    for _, count in sc:
        of.write(' {:3d}'.format(count))
    of.write('\n')

#writing the type of cell used
    if cartesian:
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

for i in range(7):
    create_POSCAR(i)
