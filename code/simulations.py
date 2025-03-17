import SierpinskiGausket.SierpinskiGausket as sg
import VicsekSet.VicsekSet as vicsekSet
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from pprint import pprint
import mpmath
import sympy
import math

np.random.seed(100)

class simulation:

    def __init__(self, n: int,  fractalConstruction, s: int):
        self.fractal = fractalConstruction
        self.n = int(n)
        self.s = int(s)


    def laplacianOperatorMatrix(self, n: int = None) -> list:
        """ 
        Finds the laplacain operators matrix. It is given by L = A - D, where

            L : is the laplacian
            A : is the adjacency matrix
            D : is the degree matrix

        Args
        ----
            n : int
                Recursion depth
        
        Returns
        -------
            list
                List of list (matrix) with degree on the diagonal and the adjacency else. It is a symmetric matrix
        """
        
        # Setup if None is given as arguments
        if n == None:
            n = self.n
        allPointsDict = self.fractal(n)
        return [[len(allPointsDict[i]["neighbours"]) if i == j else -1 if allPointsDict[j]["pos"] in allPointsDict[i]["neighbours"] else 0 for i in allPointsDict] for j in allPointsDict]

    def printLaplacianOperatorMatrix(self, n: int = None) -> None:
        """ Function to print out the laplacian in a managable way.
        
        Args
        ----
            n : int
                Recursion depth
        """

        # Setup if None is given as arguments
        if n == None:
            n = self.n

        for i in self.laplacianOperatorMatrix(n = n):
            print(i)

    def eigenVectorsAndValues(self):
        vector, matrix = mpmath.mp.eigh(mpmath.mp.matrix(self.laplacianOperatorMatrix()))
        return np.array([float(x) for x in vector]) * -(np.isclose(np.array([float(x) for x in vector], dtype=np.float64), 0) - 1), np.matrix(matrix.tolist(), dtype=np.float64) * -(np.isclose(np.matrix(matrix.tolist(), dtype=np.float64),0)-1)
    
    def npEigenVectorsAndValues(self):
        vector, matrix = np.linalg.eigh(self.laplacianOperatorMatrix())
        return vector * -(np.isclose(vector, 0) - 1), matrix * -(np.isclose(matrix, 0)-1)
    

    def DFDGsim(self, eigenfunction, s: int = None):
        n = self.n
        if s == None:
            s = self.s
        
        eigVal, eigVec = eigenfunction()
        Vn = self.fractal()

        whiteNoise = np.random.standard_normal(size = len(Vn))

        j = 0
        for point in Vn:
            Vn[point][f"X"] = sum([pow(eigVal[i],-s)*eigVec[j,i]*whiteNoise[j] if eigVal[i] > 0 else 0 for i in range(len(eigVal))])
            j +=1

        return Vn
    

    def drawThePretty(self):
        sim = self.DFDGsim(eigenfunction=self.npEigenVectorsAndValues, s = 1)
        pointx, pointy, colorValues = zip(*[[sim[x]["pos"][0], sim[x]["pos"][1], sim[x]["X"]] for x in sim])
        #pointxList, pointyList, colorValueList = list(pointx), list(pointy), list(colorValues)
        pointxList, pointyList, colorValueList = [], [], []
        
        for midpoint in [ele for ele in sim if len(sim[ele]["neighbours"]) > 2 ]:
            pointxTemp, pointyTemp, colorValuesTemp = zip(*[(np.linspace(midpoint[0], point[0]),np.linspace(midpoint[1], point[1]), np.linspace(sim[midpoint]["X"], sim[point]["X"])) for point in sim[midpoint]["neighbours"]])
            pointxList += [x for arrays in pointxTemp for x in arrays]
            pointyList += [x for arrays in pointyTemp for x in arrays]
            colorValueList += [x for arrays in colorValuesTemp for x in arrays]


        print(f"{len(pointxList)}, {len(pointyList)}, {len(colorValueList)}")
        
        tis = colorValueList + list(colorValues)
        floored = [math.floor(x) for x in tis]
        xtra = floored.sort()
        xtra2 = tis.sort()
        pprint(f"{xtra}")
        pprint(f"{tis}")

        fig, ax = plt.subplots()

        plt.scatter(pointxList, pointyList, c=cm.RdYlBu(colorValueList), s = 1)
        plt.scatter(pointx, pointy, c=cm.RdYlBu(colorValues),s = 5)

        init_length=3
        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, self.n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, self.n) * init_length - buffer, pow(3, self.n) * init_length + buffer
        y_min, y_max = -pow(3, self.n) * init_length - buffer, pow(3, self.n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        

        plt.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min(colorValues), vmax=max(colorValues)), cmap=cm.RdYlBu),ax = ax, location="right", orientation="vertical")
        plt.show()



h = 5

#sierpinski = simulation(n = h, s = 1, fractalConstruction=sg.sierpinski(n = h).pointsAndNeighbours)
#vicsek = simulation(n = h, fractalConstruction=sg.sierpinski(n=h).pointsAndNeighbours)
vicsek = simulation(n = h, s = 1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbours)

vicsek.drawThePretty()

#print(f"n=0: {vicsekSet.vicsek(n=0).pointsAndNeighbours()}")
#print(f"n=1: {vicsekSet.vicsek(n=1).pointsAndNeighbours()}")

#sim = vicsek.DFDGsim(eigenfunction=vicsek.eigenVectorsAndValues)
#im = vicsek.DFDGsim(eigenfunction=vicsek.npEigenVectorsAndValues)


#for x in sim:
#    pprint(f"{x}: {sim[x]["X"]}")
    
    #pprint(f"npsim: {npsim}")

#print(f"max: {sim[(27, 54)]["X"]}")


"""
pointx, pointy, colorValues = zip(*[(sim[x]["pos"][0], sim[x]["pos"][1], sim[x]["X"]) for x in sim])
#LineCollection(c=mpl.cm.hot())

fig, ax = mpl.pyplot.subplots()

mpl.pyplot.scatter(pointx, pointy, c=cm.viridis(colorValues/max(colorValues)), s = 8)

x = list(colorValues)
x.sort()
pprint(f"sorted: {[math.floor(y) for y in x]}")

init_length=3
n = h
# Computes limites of graph (how far x and y axis should stretch out) based on n
buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

ax.set_xlim(x_min, x_max)
ax.set_ylim(y_min, y_max)
mpl.pyplot.show()


vicsek.printLaplacianOperatorMatrix()

m, e = vicsek.eigenVectorsAndValues()
lapOp = vicsek.laplacianOperatorMatrix()

npm, npe = vicsek.npEigenVectorsAndValues()

print(f"Eigenvectors: \n {e}" )
print(f"Eigenvalues: \n {m}")

#for i in range(len(e)):
#    print(f"eigen times vec: {m[i] * e[:, i]}")
#    print(f"mat times vec: {mpmath.mp.matrix(lapOp) * e[:, i]}")


print(f"NP Eigenvectors: \n {npe}")
print(f"NP Eigenvalues: \n {npm}")
for i in range(len(npe)):
    print(f"NP eigen times vec: {npm[i] * npe[:, i]}")
    print(f"NP mat times vec: {lapOp @ npe[:, i]}")

for i in npe:
    print(i)

#print(mpmath.mp.inverse(e))
#print(e)


#print(0.8660254 * k)"
"""