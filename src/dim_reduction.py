import MDAnalysis as mda
from MDAnalysis.analysis import pca
import nglview as nv
import matplotlib.pylab as plt
import pandas as pd
import seaborn as sns
#constants 
PDB_ROOT = "./data/1ULL.pdb"
FILE_NAME = "dim_reduction"

u = mda.Universe("./data/1ULL.pdb", dt=1.0)
# print(u.residues) #just to make sure it's readign

pc = pca.PCA(u, select="backbone", align=True, mean=None, n_components=None)
pc.run()

backbone = u.select_atoms('backbone')
n_bb = len(backbone)
print(f"Number of backbone atoms:{n_bb}")
print(f"Shape: {pc.results.p_components.shape}") 
print(f"Cumulative Variance %: {pc.results.cumulated_variance}")

#plot
plt.plot(pc.cumulated_variance[:10])
plt.title("Line Chart of Cumulated Variance")
plt.xlabel('Principal component')
plt.ylabel('Cumulative variance')
# plt.show()
plt.savefig(f"results/PCA_cumulative_variance.png")

#transform
N_COMPONENTS = 3
transformed = pc.transform(backbone, n_components=N_COMPONENTS)
# print(transformed.shape)

df = pd.DataFrame(transformed,
                  columns=['PC{}'.format(i+1) for i in range(N_COMPONENTS)])
df['Time (ps)'] = df.index * u.trajectory.dt

g = sns.PairGrid(df, hue='Time (ps)', palette=sns.color_palette('Oranges_d', n_colors=len(df)))
g.map(plt.scatter, marker='.')
# plt.show()
plt.savefig(f"results/PCA_coords.png")
