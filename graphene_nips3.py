import supercell_core as sc
import matplotlib.pyplot as plt
import numpy as np

# Read graphene and NiPS3 definition from POSCAR
#graphene = sc.read_POSCAR("/global/cscratch1/sd/debajitc/supercell-core/supercell_core/resources/vasp/graphene/POSCAR_with_atomic_species_names")
#nips3 = sc.read_POSCAR("/global/cscratch1/sd/debajitc/supercell-core/supercell_core/resources/vasp/NiPS3/POSCAR", atomic_species=['Ni', 'P', 'S'])
graphene = sc.read_POSCAR("global/cscratch1/sd/debajitc/calculations/tmds/mos2/single-layer/DFT-D3/CONTCAR")
nips3 = sc.read_POSCAR("global/cscratch1/sd/debajitc/calculations/tmds/mos2/single-layer/DFT-D3/CONTCAR", atomic_species=['Mo', 'S'])
h = sc.heterostructure().set_substrate(graphene)\
.add_layer(nips3)

#res = h.opt(max_el=2)
# thetas=\
#[np.arange(0, 7*sc.DEGREE, 0.1*sc.DEGREE)])

#res.superlattice().save_POSCAR("POSCAR_sc")

res = h.opt(max_el=8, thetas=[
np.arange(0, 30 * sc.DEGREE, 0.25*sc.DEGREE)])

res.superlattice().save_POSCAR("POSCAR_sc")

# Draw the resulting supercell
res.superlattice().draw()
plt.title("""Best supercell found for $\\theta$ = 
{} max strain = {:.4g}"""
.format(res.thetas()[0], res.max_strain()))
