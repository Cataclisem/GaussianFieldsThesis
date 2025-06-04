import SierpinskiGausket.SierpinskiGausket as sg
import VicsekSet.VicsekSet as vicsekSet
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from pprint import pprint
import mpmath
from scipy import sparse
from scipy.sparse import linalg
import math
import timeit

import multiprocessing as mp

#np.random.seed(100)

class simulation:

    def __init__(self, n: int,  fractalConstruction, s: int, fractalType: str):
        self.fractal = fractalConstruction
        self.n = int(n)
        self.s = int(s)
        self.SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

        if fractalType.lower() in {"sierpinski", "vicsek"}:
            self.fractalType = fractalType.lower()
            self.pointAmount = (pow(5, self.n) *4 +1)*(self.fractalType == "vicsek") + (int(3*(pow(3,self.n) + 1)/2))*(self.fractalType=="sierpinski")
        else:
            raise Exception("Did you mean 'vicsek' or 'sierpinski'")
    
    def pointAmont(self):
        return self.pointAmont

    def laplacianOperatorMatrix(self, n: int = None) -> list:
        """ 
        Finds the laplacain operators matrix. It is given by L = A - D, where

            L : is the laplacian
            A : is the adjacency matrix
            D : is the degree matrix
        
        This function is now multiprocessing, since this has to go fast.
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
        print(f"Start Laplacian {int(timeStart)}")
        allPointsDict = self.fractal(n)
        pool = mp.Pool(mp.cpu_count()) #Initialize the multiprocessing
        who = pool.starmap(self.laplacianOperatorMatrixFunction, [([x], allPointsDict) for x in allPointsDict])
        timeEnd = timeit.default_timer()
        print(f"end Laplacian {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        return who
    
    
    def laplacianOperatorMatrixFunction(self, i, allPointsDict: dict):
        """ This is a function for multiprocessing
        
        """
        return [len(allPointsDict[i[0]]["neighbours"]) if i[0] == j else -1 if j in allPointsDict[i[0]]["neighbours"] else 0 for j in allPointsDict]



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
        print(f"Start Eigen {int(timeStart)}")
        vector, matrix = np.linalg.eigh(self.laplacianOperatorMatrix())
        timeEnd = timeit.default_timer()
        print(f"end Eigen {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        return vector * -(np.isclose(vector, 0) - 1), matrix * -(np.isclose(matrix, 0)-1)
    
    def scipyEigenVectorsAndValues(self):
        """ scipys sparse, cant find all eigenvectors but it can find close to all. We can therefore use the fact that the graph is connected to know that there is 1 zero eigenvalue. By some miracle it is also the one scipy doesn't find.
        """
        timeStart =timeit.default_timer()
        print(f"Start Scipy Eigen {int(timeStart)}")
        vector, matrix = linalg.eigsh(sparse.csc_matrix(self.laplacianOperatorMatrix()), k=self.pointAmount - 1)
        timeEnd = timeit.default_timer()
        print(f"end scipy Eigen {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        return vector, matrix
    

    def DFDGsim(self, eigenfunction, s: int = None, whiteNoise: np.ndarray = None):
        if s == None:
            s = self.s
        
        eigVal, eigVec = eigenfunction()
        print(f"eigenVec: {eigVec}")
        print(f"eigenvec: {eigVec[0,1]}")
        Vn = self.fractal()

        if not isinstance(whiteNoise, np.ndarray):
            whiteNoise = np.random.standard_normal(size = len(Vn))

        j = 0
        for point in Vn:
            Vn[point][f"X"] = sum([pow(eigVal[i],-s)*eigVec[j,i]*whiteNoise[i] if eigVal[i] > 0 else 0 for i in range(len(eigVal))])
            j +=1
            if j % 1000 == 0:
                print(f"Point reached: {j}")

        return Vn
    

    def MakeThePretty(self, s: int, whiteNoise: list = None):
        sim = self.DFDGsim(eigenfunction=self.npEigenVectorsAndValues, s = s, whiteNoise=whiteNoise)
        pointx, pointy, colorValues = zip(*[[x[0], x[1], sim[x]["X"]] for x in sim])
        pointxList, pointyList, colorValuesList = [], [], []
        
        if self.fractalType == "vicsek":
            for midpoint in [ele for ele in sim if len(sim[ele]["neighbours"]) > 2 ]:
                pointxTemp, pointyTemp, colorValuesTemp = zip(*[(np.linspace(midpoint[0], point[0]),np.linspace(midpoint[1], point[1]), np.linspace(sim[midpoint]["X"], sim[point]["X"])) for point in sim[midpoint]["neighbours"]])
                pointxList += [x for arrays in pointxTemp for x in arrays]
                pointyList += [x for arrays in pointyTemp for x in arrays]
                colorValuesList += [x for arrays in colorValuesTemp for x in arrays]

        elif self.fractalType == "sierpinski":
            checkedPoints = set()
      
            for midpoint in sim:
                for point in [x for x in sim[midpoint]["neighbours"] if (midpoint, x) not in checkedPoints]:
                    pointxTemp, pointyTemp, colorValuesTemp = zip(*[(np.linspace(midpoint[0], point[0]),np.linspace(midpoint[1], point[1]), np.linspace(sim[midpoint]["X"], sim[point]["X"])) for point in sim[midpoint]["neighbours"]])
                    pointxList += [x for arrays in pointxTemp for x in arrays]
                    pointyList += [x for arrays in pointyTemp for x in arrays]
                    colorValuesList += [x for arrays in colorValuesTemp for x in arrays]
                    checkedPoints.update({(midpoint, point)}, {(point, midpoint)})

        print(f"Total amount of points: {len(pointxList)}, {len(pointyList)}, {len(colorValuesList)} \n")

        return pointx, pointy, colorValues, pointxList, pointyList, colorValuesList
        
    
    def drawThePretty(self, sValues: list, sameWhiteNoise: bool = False, nrows: int = 2, whiteNoise = None):
        
        if sameWhiteNoise == True:
            whiteNoise = np.random.standard_normal(size = self.pointAmount)
        elif whiteNoise.any() != None:
            whiteNoise = whiteNoise
        else:
            whiteNoise = None
        
        print(f"WhiteNoise: {whiteNoise}")

        if len(sValues) > 2:
            sValHalfRdUp = math.ceil(len(sValues)/nrows) 
            fig, axes = plt.subplots(nrows= nrows, ncols=sValHalfRdUp)
            for i in range(nrows):
                for ax, s in zip(axes[i], sValues[i*sValHalfRdUp:i*sValHalfRdUp + sValHalfRdUp]):
                    self.genGraph(fig = fig, ax=ax, s=s, whiteNoise=whiteNoise)
        elif len(sValues) == 2:
            fig, axes = plt.subplots(ncols=2)
            for ax, s in zip(axes, sValues):
                self.genGraph(fig = fig, ax=ax, s=s, whiteNoise=whiteNoise)
        else:
            fig, ax = plt.subplots()
            self.genGraph(fig=fig, ax=ax, s=sValues[0], whiteNoise=whiteNoise)


    def genGraph(self, fig: plt.Figure, ax: plt.Axes, s: int, whiteNoise: np.ndarray):
        print(f"Starting on s={s}")
        pointx, pointy, colorValues, pointxList, pointyList, colorValuesList = self.MakeThePretty(s=s, whiteNoise=whiteNoise)
        ax.set_title(f"V{self.n}: s = {s}")
        ax.scatter(pointxList, pointyList, c=cm.brg(colorValuesList), s = 1)
        ax.scatter(pointx, pointy, c=cm.brg(colorValues), s = 3)
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
        ax.set_axis_off()
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min(colorValues), vmax=max(colorValues)), cmap=cm.brg),ax = ax, location="right", orientation="vertical")



h = 4

sierpinski = simulation(n = h, s = 1, fractalConstruction=sg.sierpinski(n = h).pointsAndNeighbours, fractalType="sierpinski")
#vicsek = simulation(n = h, fractalConstruction=sg.sierpinski(n=h).pointsAndNeighbours)
vicsek = simulation(n = h, s = 1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbourswhat, fractalType="vicsek")

#needs = sg.sierpinski(n=h).pointsAndNeighbours()
#for x in needs:
#   print(f"{x}: {needs[x]}")


if __name__ == '__main__':    
    print(f"h : {h}")
    #vicsek.printLaplacianOperatorMatrix()
    #vec, mat = vicsek.scipyEigenVectorsAndValues()
    #npm, npe = vicsek.npEigenVectorsAndValues()

    #print(f"Scipy vec: \n {vec}")
    #print(f"NP Eigenvalues: \n {npm}")

    #for i in range(len(vec)):
    #    if not np.isclose(vec[i], npm[i+1]) or vec[i] == 0:
    #        print(f"sci: {vec[i]}, NP: {npm[i+1]}")

    #for x in mat:
    #    print(x)

    whatToRun = 3

    dh = math.log(5)/math.log(3)
    dw = dh +1

    if whatToRun == 1:
        print(f"# points: {sierpinski.pointAmount}")
        start = timeit.default_timer()
        print("Sier")
        sierpinski.drawThePretty(sValues=[1], sameWhiteNoise=False)
        end = timeit.default_timer()
        print(f"It took {end - start} seconds")

    if whatToRun == 2:
        print(f"Starting process with: \n # Vicsek points: {vicsek.pointAmount}")
        start = timeit.default_timer()
        vicsek.drawThePretty(sValues=[(2* (math.log(5)/math.log(3)))/(1+(math.log(5)/math.log(3))), (2* (math.log(5)/math.log(3)))/(1+(math.log(5)/math.log(3)))])
        vicsek.drawThePretty(sValues=[(math.log(5)/math.log(3))/(2* (1+(math.log(5)/math.log(3)))), (math.log(5)/math.log(3))/(2* (1+(math.log(5)/math.log(3))))])
        #vicsek.drawThePretty(sValues=[0.001, 0.001], sameWhiteNoise=False)
        #vicsek.drawThePretty(sValues=[1,1], sameWhiteNoise=False)
        #vicsek.drawThePretty(sValues=[5,5], sameWhiteNoise=False)
        #vicsek.drawThePretty(sValues=[20,20], sameWhiteNoise=False)
        end = timeit.default_timer()
        print(f"It took {end - start} seconds")
    
    if whatToRun == 3:
        i=0
        for h in [5, 5]: 
            WN = np.random.standard_normal(pow(5, h)*4 + 1)
            for sVals in [0.001, 0.5, 1, 20]:
                print(f"Vm: {h}")
                #theOne = simulation(n = h, s = dh/(2*dw) +0.1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbourswhat, fractalType="vicsek")
                theOne = simulation(n = h, s = dh/(2*dw) +0.1, fractalConstruction=sg.sierpinski(n=h).pointsAndNeighbours, fractalType="sierpinski")
                theOne.drawThePretty(sValues=[sVals], whiteNoise=WN)
                #theOne.drawThePretty(sValues=[0.001, 0.5, 1, 20], sameWhiteNoise=True, whiteNoise=)
                #plt.gca().set_position([0, 0, 1, 1])
                plt.savefig(f"c:/Users/chris/GaussianFieldsThesis/code/sim/sierpinski/sameWhite/{theOne.fractalType}Sim_V{h}_SameWhiteNoise{i}-s{str(sVals).replace(".","_")}.svg", format="svg")
                i+=1
    #plt.show()
