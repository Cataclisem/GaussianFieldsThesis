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
import timeit

np.random.seed(100)

class simulation:

    def __init__(self, n: int,  fractalConstruction, s: int, fractalType: str):
        self.fractal = fractalConstruction
        self.n = int(n)
        self.s = int(s)

        if fractalType.lower() in {"sierpinski", "vicsek"}:
            self.fractalType = fractalType.lower()
            self.pointAmount = (pow(5, self.n) *4 +1)*(self.fractalType == "vicsek") + (int(3*(pow(3,self.n) + 1)/2))*(self.fractalType=="sierpinski")
        else:
            raise Exception("Did you mean 'vicsek' or 'sierpinski'")


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
        timeStart =timeit.default_timer()
        print(f"Start Laplacian {timeStart}")
        allPointsDict = self.fractal(n)
        timeEnd = timeit.default_timer()
        print(f"end Laplacian {timeEnd}, and i took: {timeEnd - timeStart}")
        return [[len(allPointsDict[i]["neighbours"]) if i == j else -1 if j in allPointsDict[i]["neighbours"] else 0 for i in allPointsDict] for j in allPointsDict]

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
        timeStart =timeit.default_timer()
        print(f"Start Eigen {timeStart}")
        vector, matrix = np.linalg.eigh(self.laplacianOperatorMatrix())
        timeEnd = timeit.default_timer()
        print(f"end Eigen {timeEnd}, and i took: {timeEnd - timeStart}")
        return vector * -(np.isclose(vector, 0) - 1), matrix * -(np.isclose(matrix, 0)-1)
    

    def DFDGsim(self, eigenfunction, s: int = None, whiteNoise: np.ndarray = None):
        n = self.n
        if s == None:
            s = self.s
        
        eigVal, eigVec = eigenfunction()
        Vn = self.fractal()

        if not isinstance(whiteNoise, np.ndarray):
            whiteNoise = np.random.standard_normal(size = len(Vn))

        j = 0
        for point in Vn:
            Vn[point][f"X"] = sum([pow(eigVal[i],-s)*eigVec[j,i]*whiteNoise[j] if eigVal[i] > 0 else 0 for i in range(len(eigVal))])
            j +=1
            if j % 1000 == 0:
                print(f"j: {j}")

        #for x in Vn:
        #    print(f"allDict: {x}: {Vn[x]}")

        return Vn
    

    def MakeThePretty(self, s: int, whiteNoise: list = None):
        sim = self.DFDGsim(eigenfunction=self.npEigenVectorsAndValues, s = s, whiteNoise=whiteNoise)
        pointx, pointy, colorValues = zip(*[[x[0], x[1], sim[x]["X"]] for x in sim])
        #pointxList, pointyList, colorValueList = list(pointx), list(pointy), list(colorValues)
        pointxList, pointyList, colorValuesList = [], [], []
        
        if self.fractalType == "vicsek":
            for midpoint in [ele for ele in sim if len(sim[ele]["neighbours"]) > 2 ]:
                pointxTemp, pointyTemp, colorValuesTemp = zip(*[(np.linspace(midpoint[0], point[0]),np.linspace(midpoint[1], point[1]), np.linspace(sim[midpoint]["X"], sim[point]["X"])) for point in sim[midpoint]["neighbours"]])
                pointxList += [x for arrays in pointxTemp for x in arrays]
                pointyList += [x for arrays in pointyTemp for x in arrays]
                colorValuesList += [x for arrays in colorValuesTemp for x in arrays]

        elif self.fractalType == "sierpinski":
            checkedPoints = set()
            #for x in sim:
                #pprint(f"{x}: {sim[x]}")
            for midpoint in sim:
                for point in [x for x in sim[midpoint]["neighbours"] if (midpoint, x) not in checkedPoints]:
                    pointxTemp, pointyTemp, colorValuesTemp = zip(*[(np.linspace(midpoint[0], point[0]),np.linspace(midpoint[1], point[1]), np.linspace(sim[midpoint]["X"], sim[point]["X"])) for point in sim[midpoint]["neighbours"]])
                    pointxList += [x for arrays in pointxTemp for x in arrays]
                    pointyList += [x for arrays in pointyTemp for x in arrays]
                    colorValuesList += [x for arrays in colorValuesTemp for x in arrays]
                    checkedPoints.update({(midpoint, point)}, {(point, midpoint)})
                    #print(f"midpoint: {midpoint} \n Neigh point: {point} \n checked: {checkedPoints}")



        print(f"{len(pointxList)}, {len(pointyList)}, {len(colorValuesList)}")
        #print(f"{s}: {colorValues}")

        return pointx, pointy, colorValues, pointxList, pointyList, colorValuesList
        
    
    def drawThePretty(self, sValues: list, sameWhiteNoise: bool = False, nrows: int = 2):
        
        if sameWhiteNoise == True:
            whiteNoise = np.random.standard_normal(size = self.pointAmount)
        else:
            whiteNoise = None

        if len(sValues) > 1:
            sValHalfRdUp = math.ceil(len(sValues)/nrows) 
            fig, axes = plt.subplots(nrows= nrows, ncols=sValHalfRdUp)
            for i in range(nrows):
                for ax, s in zip(axes[i], sValues[i*sValHalfRdUp:i*sValHalfRdUp + sValHalfRdUp]):
                    self.genGraph(fig = fig, ax=ax, s=s, whiteNoise=whiteNoise)
        else:
            fig, ax = plt.subplots()
            self.genGraph(fig=fig, ax=ax, s=sValues[0], whiteNoise=whiteNoise)


    def genGraph(self, fig: plt.Figure, ax: plt.Axes, s: int, whiteNoise: np.ndarray):
        pointx, pointy, colorValues, pointxList, pointyList, colorValuesList = self.MakeThePretty(s=s, whiteNoise=whiteNoise)
        ax.set_title(f"s = {s}")
        ax.scatter(pointxList, pointyList, c=cm.brg(colorValuesList), s = 1)
        ax.scatter(pointx, pointy, c=cm.brg(colorValues), s = 20)
        # Computes limites of graph (how far x and y axis should stretch out) based on n
        if self.fractalType == "vicsek":
            init_length=3
            buffer = (pow(3, self.n) * init_length)/10 # A buffer to make the final graph not look as cramped
            x_min, x_max = -pow(3, self.n) * init_length - buffer, pow(3, self.n) * init_length + buffer
            y_min, y_max = -pow(3, self.n) * init_length - buffer, pow(3, self.n) * init_length + buffer
        elif self.fractalType == "sierpinski":
            buffer = pow(2, self.n)/10 # A buffer to make the final graph not look as cramped
            x_min, x_max = -buffer, pow(2, self.n) + buffer
            y_min, y_max = -buffer, pow(2, self.n) + buffer
        
        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min(colorValues), vmax=max(colorValues)), cmap=cm.brg),ax = ax, location="right", orientation="vertical")



h = 7

sierpinski = simulation(n = h, s = 1, fractalConstruction=sg.sierpinski(n = h).pointsAndNeighbours, fractalType="sierpinski")
#vicsek = simulation(n = h, fractalConstruction=sg.sierpinski(n=h).pointsAndNeighbours)
vicsek = simulation(n = h, s = 1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbourswhat, fractalType="vicsek")

#needs = sg.sierpinski(n=h).pointsAndNeighbours()
#for x in needs:
#   print(f"{x}: {needs[x]}")

whatToRun = 2

if whatToRun == 1:
    print(f"# points: {sierpinski.pointAmount}")
    start = timeit.default_timer()
    print("Sier")
    sierpinski.drawThePretty(sValues=[1], sameWhiteNoise=True)
    end = timeit.default_timer()
    print(f"It took {end - start} seconds")

if whatToRun == 2:
    print(f"# points: {vicsek.pointAmount}")
    start = timeit.default_timer()
    print("vic")
    vicsek.drawThePretty(sValues=[1], sameWhiteNoise=True)
    end = timeit.default_timer()
    print(f"It took {end - start} seconds")
plt.show()

#vicsek.drawThePrettyVicsek(sValues = [0.01, 1, 2, 10], sameWhiteNoise=True)

#plt.show()
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