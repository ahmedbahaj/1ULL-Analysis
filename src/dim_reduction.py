import MDAnalysis as mda
from MDAnalysis.analysis import pca
import nglview as nv

#constants 
PDB_ROOT = "./data/1ULL.pdb"


u = mda.Universe("./data/1ULL.pdb")
# print(u.residues) #just to make sure it's readign

pc = pca.PCA(u, select="backbone", align=True, mean=None, n_components=None)
pc.run()

backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print("number of backbone :{n_bb}")
print("shape: {pc.p_components.shape}") 