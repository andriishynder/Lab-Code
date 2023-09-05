# %%
import math
import numpy as np
from numpy.linalg import eig
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# %%
## Formulas and definitions

def MassCentre(molecule):
    moleculeMass = 0
    for m in molecule.values():
        moleculeMass = moleculeMass + m["mass"]
    
    resCentre = [0, 0, 0]
    for m in molecule.values():
        resCentre[0] = resCentre[0] + m["mass"] * m["crd"][0]
        resCentre[1] = resCentre[1] + m["mass"] * m["crd"][1]
        resCentre[2] = resCentre[2] + m["mass"] * m["crd"][2]

    resCentre[0] = round(resCentre[0] / moleculeMass, 4)
    resCentre[1] = round(resCentre[1] / moleculeMass, 4)
    resCentre[2] = round(resCentre[2] / moleculeMass, 4)    

    return resCentre

def KroneckerDeltaFunc(i, j):
    if(i == j): return 1
    else: return 0

def I_func(i, j, molecule, centre):
    res = 0
    for m in molecule.values():
        res = res + m["mass"] * (math.dist(centre, m["crd"])**2 * KroneckerDeltaFunc(i, j) - (m["crd"][i] - centre[i])*(m["crd"][j] - centre[i]))
    return res    

def InertiaTensor(molecule, centre):
    resI = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            resI[i][j] = I_func(i,j, molecule, centre)
    return resI

def ChangeBasis(molecule, newOrigin):
        for m in molecule.values():
            for i in range(3):
                m["crd"][i] = m["crd"][i] - newOrigin[i]

# %%

## Molecule definition

atom_C = { "crd": [0.0015, -0.7180, -0.0350], "mass":12.0}
atom_O = { "crd": [-0.0019, 0.6780, 0.0835], "mass":15.995}
atom_H0 = { "crd": [0.0058, 1.0768, -0.7672], "mass":1.0078}
atom_H1 = { "crd": [-0.8762, -1.0812, -0.5626], "mass":1.0078}
atom_H2 = { "crd": [-0.0081, -1.1248, 0.9656], "mass":1.0078}
atom_H3 = { "crd": [0.8909, -1.0784, -0.5446], "mass":1.0078}
molecule = {
    "C": atom_C, "O": atom_O,  "H0": atom_H0, "H1": atom_H1, "H2": atom_H2, "H3": atom_H3
    }

moleculaPos = list(map(lambda m: m["crd"], molecule.values()))

moleculeMassCentre = MassCentre(molecule)
print("Centre of mass : {}".format(moleculeMassCentre)) 

ChangeBasis(molecule, moleculeMassCentre)
moleculaPos = list(map(lambda m: m["crd"], molecule.values()))
print("Coordinates in 'centre of mass' basis: {}".format(moleculaPos))

moleculeMassCentre = MassCentre(molecule)
#print("Centre of mass : {}".format(moleculeMassCentre)) 

# %%
I = InertiaTensor(molecule, [0, 0, 0])

print("Tensor of intertia:") 
#print(np.matrix(I))
print(I)

# %%

# Eigen - decomposition
evalues,evectors = eig(I)
#print(evalues)
#print(evectors)

D = np.zeros((3,3))
for i in range(3):
    D[i,i] = evalues[i]

print("Prinicpal moments of intertia: {}".format(D))

u1 = [row[0] for row in evectors]
u2 = [row[1] for row in evectors]
u3 = [row[2] for row in evectors]

print("Principal axes:")
print("u1 : {}".format(u1))
print("u2 : {}".format(u2))
print("u3 : {}".format(u3))

# %%

## Drawing definition

def DrawMolecule(molecule, animate=False):
    def link(a1,a2,i):
        return [a1["crd"][i], a2["crd"][i]]
    
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    # Draw atoms
    ax.scatter([molecule["C"]["crd"][0]],[molecule["C"]["crd"][1]],[molecule["C"]["crd"][2]], c="green", s=700)
    ax.scatter([molecule["O"]["crd"][0]],[molecule["O"]["crd"][1]],[molecule["O"]["crd"][2]], c="red", s=700)
    ax.scatter([molecule["H0"]["crd"][0]],[molecule["H0"]["crd"][1]],[molecule["H0"]["crd"][2]], c="black", s=300, marker="o")
    ax.scatter([molecule["H1"]["crd"][0]],[molecule["H1"]["crd"][1]],[molecule["H1"]["crd"][2]], c="black", s=300, marker="o")
    ax.scatter([molecule["H2"]["crd"][0]],[molecule["H2"]["crd"][1]],[molecule["H2"]["crd"][2]], c="black", s=300, marker="o")
    ax.scatter([molecule["H3"]["crd"][0]],[molecule["H3"]["crd"][1]],[molecule["H3"]["crd"][2]], c="black", s=300, marker="o")

    # Draw links
    linkWidth = 10
    ax.plot(link(molecule["C"], molecule["O"], 0),link(molecule["C"], molecule["O"], 1),link(molecule["C"], molecule["O"], 2), c='darkgray', linewidth = linkWidth)
    ax.plot(link(molecule["O"], molecule["H0"], 0),link(molecule["O"], molecule["H0"], 1),link(molecule["O"], molecule["H0"], 2), c='darkgray', linewidth = linkWidth)
    ax.plot(link(molecule["C"], molecule["H1"], 0),link(molecule["C"], molecule["H1"], 1),link(molecule["C"], molecule["H1"], 2), c='darkgray', linewidth = linkWidth)
    ax.plot(link(molecule["C"], molecule["H2"], 0),link(molecule["C"], molecule["H2"], 1),link(molecule["C"], molecule["H2"], 2), c='darkgray', linewidth = linkWidth)
    ax.plot(link(molecule["C"], molecule["H3"], 0),link(molecule["C"], molecule["H3"], 1),link(molecule["C"], molecule["H3"], 2), c='darkgray', linewidth = linkWidth)

    # Draw mass centre
    moleculeMassCentre = MassCentre(molecule)
    ax.scatter([moleculeMassCentre[0]],[moleculeMassCentre[1]],[moleculeMassCentre[2]], c="blue", s=30, marker="^")

    # Draw coordinate axis
    ax.scatter([0],[0],[0], c='black')

    ax.plot([0, 1], [0, 0], [0, 0], c='black', linewidth = 0.5)
    ax.plot([0, 0], [0, 1], [0, 0], c='black', linewidth = 0.5)
    ax.plot([0, 0], [0, 0], [0, 1], c='black', linewidth = 0.5)

    # Draw eigenvector basis
    ax.plot([0, u1[0]], [0, u1[1]], [0, u1[2]], c='red', linewidth = 0.5)
    ax.text(u1[0], u1[1], u1[2], "x'", c='red')
    ax.plot([0, u2[0]], [0, u2[1]], [0, u2[2]], c='red', linewidth = 0.5)
    ax.text(u2[0], u2[1], u2[2], "y'", c='red')
    ax.plot([0, u3[0]], [0, u3[1]], [0, u3[2]], c='red', linewidth = 0.5)
    ax.text(u3[0], u3[1], u3[2], "z'", c='red')

    # Animate rotation
    if animate:
        for angle in range(0, 360):
            ax.view_init(30, angle)
            plt.draw()
            plt.pause(.001)
    else:
        plt.show()

# %%
DrawMolecule(molecule)