import sys
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

    def __init__(self, n: int,  fractalConstruction, s: float, fractalType: str):
        self.fractal = fractalConstruction
        self.n = int(n)
        self.s = float(s)
        self.SUB = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")

        if fractalType.lower() in {"sierpinski", "vicsek"}:
            self.fractalType = fractalType.lower()
            self.pointAmount = (pow(5, self.n) *4 +1)*(self.fractalType == "vicsek") + (int(3*(pow(3,self.n) + 1)/2))*(self.fractalType=="sierpinski")
        else:
            raise Exception("Did you mean 'vicsek' or 'sierpinski'")
    
    def pointAmont(self):
        return self.pointAmont

    def laplacianOperatorMatrix(self, n: int = None, withConst: bool = False) -> list:
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
        point_order = list(allPointsDict.keys())
        pool = mp.Pool(mp.cpu_count()) #Initialize the multiprocessing
        if withConst == True:
            who = pool.starmap(self.laplacianOperatorMatrixFunctionWithConstants, [([x], allPointsDict) for x in allPointsDict])
        else:
            who = pool.starmap(self.laplacianOperatorMatrixFunction, [([x], allPointsDict) for x in allPointsDict])
        timeEnd = timeit.default_timer()
        print(f"end Laplacian {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        return who
    
    
    def laplacianOperatorMatrixFunctionWithConstants(self, i, allPointsDict: dict):
        """ This is a function for multiprocessing
        
        """
        if self.fractalType == "vicsek":
            print("IN VICSEK LAPLACIAN")
            if len(allPointsDict[i[0]]["neighbours"]) == 1:
                const = (3**self.n) *(pow(5, -self.n - 1)**(-1))
            elif len(allPointsDict[i[0]]["neighbours"]) == 2:
                const = (3**self.n) *((2*pow(5, -self.n - 1))**(-1))
            elif len(allPointsDict[i[0]]["neighbours"]) == 4:
                const = (3**self.n) * (pow(5, -self.n)**(-1))
            else:
                const = 1
        
        elif self.fractalType == "sierpinski":
            print("IN SIERPINSKI LAPLACIAN")
            if len(allPointsDict[i[0]]["neighbours"]) == 2:
                const = 3*pow(5, self.n)
            elif len(allPointsDict[i[0]]["neighbours"]) == 4:
                const = (3/2)*pow(5, self.n)
            else:
                const = 1
        else:
            print("IN NO LAPLACIAN")
            const = 1

        return [-const*len(allPointsDict[i[0]]["neighbours"]) if i[0] == j else 1*const if j in allPointsDict[i[0]]["neighbours"] else 0 for j in allPointsDict]


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
    
    def npEigenVectorsAndValues(self, withConst: bool = False):
        timeStart =timeit.default_timer()
        print(f"Start Eigen {int(timeStart)}")
        vector, matrix = np.linalg.eigh(self.laplacianOperatorMatrix(withConst = withConst))
        timeEnd = timeit.default_timer()
        print(f"end Eigen {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        #pprint(f"vec: {vector}, isclose: {vector[0]*-(np.isclose(vector[0], 0)-1)}")
        return vector * -(np.isclose(vector, 0) - 1), matrix * -(np.isclose(matrix, 0)-1)
    

    def partwiseNPEigenVecAndValues(self):
        timeStart =timeit.default_timer()
        print(f"Start Eigen Partwise: {int(timeStart)}")
        laplacian = self.laplacianOperatorMatrix(withConst=True)
        vector, matrix = np.linalg.eigh(laplacian)
        timeEnd = timeit.default_timer()
        print(f"end Eigen {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        #pprint(f"vec: {vector}, isclose: {vector[0]*-(np.isclose(vector[0], 0)-1)}")
        return laplacian, vector, matrix, vector * -(np.isclose(vector, 0) - 1), matrix * -(np.isclose(matrix, 0)-1), self.fractal(self.n)
    
    def scipyEigenVectorsAndValues(self):
        """ scipys sparse, cant find all eigenvectors but it can find close to all. We can therefore use the fact that the graph is connected to know that there is 1 zero eigenvalue. By some miracle it is also the one scipy doesn't find.
        """
        timeStart =timeit.default_timer()
        print(f"Start Scipy Eigen {int(timeStart)}")
        vector, matrix = linalg.eigsh(sparse.csc_matrix(self.laplacianOperatorMatrix()), k=self.pointAmount - 1)
        timeEnd = timeit.default_timer()
        print(f"end scipy Eigen {int(timeEnd)}, and it took: {int(timeEnd - timeStart)} seconds")
        return vector, matrix
    

    def DFDGsim(self, eigenfunction, s: int = None, whiteNoise: np.ndarray = None, withConst: bool = False):
        if s == None:
            s = self.s
        
        eigVal, eigVec = eigenfunction(withConst = withConst)
        Vn = self.fractal()

        if not isinstance(whiteNoise, np.ndarray):
            whiteNoise = np.random.standard_normal(size = len(Vn))

        j = 0
        for point in Vn:
            print(f"-s: {-s}, point: {point}")
            for i in range(len(eigVal)):
                if eigVal[i] > 0:
                    print(f"eig: {pow(eigVal[i],-s)}, vec: {eigVec[j,i]}, WN: {whiteNoise[i]}, sum: {sum([pow(eigVal[k],-s)*eigVec[j,k]*whiteNoise[k] if eigVal[k] > 0 else 0 for k in range(len(eigVal))])}")

            Vn[point][f"X"] = sum([pow(eigVal[i],-s)*eigVec[j,i]*whiteNoise[i] if eigVal[i] > 0 else 0 for i in range(len(eigVal))])
            j +=1
            if j % 1000 == 0:
                print(f"Point reached: {j}")

        #print(f"Vn: {Vn}")
        return Vn
    

    def MakeThePretty(self, s: int, whiteNoise: list = None, withConst: bool = False):
        sim = self.DFDGsim(eigenfunction=self.npEigenVectorsAndValues, s = s, whiteNoise=whiteNoise, withConst=withConst)
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
        
    
    def drawThePretty(self, sValues: list, sameWhiteNoise: bool = False, nrows: int = 2, whiteNoise: np.ndarray = None, withConst: bool = False):
        
        print(f"isin: {isinstance(whiteNoise, np.ndarray)}")
        if sameWhiteNoise == True:
            whiteNoise = np.random.standard_normal(size = self.pointAmount)
        elif isinstance(whiteNoise, np.ndarray):
            whiteNoise = whiteNoise
        else:
            whiteNoise = None
        
        #print(f"WhiteNoise: {whiteNoise}")

        if len(sValues) > 2:
            sValHalfRdUp = math.ceil(len(sValues)/nrows) 
            fig, axes = plt.subplots(nrows= nrows, ncols=sValHalfRdUp)
            for i in range(nrows):
                for ax, s in zip(axes[i], sValues[i*sValHalfRdUp:i*sValHalfRdUp + sValHalfRdUp]):
                    self.genGraph(fig = fig, ax=ax, s=s, whiteNoise=whiteNoise, withConst = withConst)
        elif len(sValues) == 2:
            fig, axes = plt.subplots(ncols=2)
            for ax, s in zip(axes, sValues):
                self.genGraph(fig = fig, ax=ax, s=s, whiteNoise=whiteNoise, withConst = withConst)
        else:
            fig, ax = plt.subplots()
            self.genGraph(fig=fig, ax=ax, s=sValues[0], whiteNoise=whiteNoise, withConst = withConst)


    def genGraph(self, fig: plt.Figure, ax: plt.Axes, s: int, whiteNoise: np.ndarray, withConst: bool = False):
        print(f"Starting on s={s}")
        pointx, pointy, colorValues, pointxList, pointyList, colorValuesList = self.MakeThePretty(s=s, whiteNoise=whiteNoise, withConst = withConst)
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
        ax.set_ylim(y_min, y_max) # type: ignore
        ax.set_axis_off()
        fig.colorbar(cm.ScalarMappable(norm=colors.Normalize(vmin=min(colorValues), vmax=max(colorValues)), cmap=cm.brg),ax = ax, location="right", orientation="vertical")



h = 1

sierpinski = simulation(n = h, s = 1, fractalConstruction=sg.sierpinski(n = h).pointsAndNeighbours, fractalType="sierpinski")
vicsek = simulation(n = h, s = 1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbourswhat, fractalType="vicsek")

if __name__ == '__main__':    
    print(f"h : {h}")
    whatToRun = 1

    dh = math.log(5)/math.log(3)
    dw = dh +1
    
    def printMatrix(matrix):
        for i in matrix:
            print(i)

    if whatToRun == 1:
        np.set_printoptions(threshold=sys.maxsize)
        White_Noise = np.random.standard_normal(pow(5, h)*4 + 1)
        laplacian, OGvec, OGmat, vec, mat, allPoints = vicsek.partwiseNPEigenVecAndValues()
        #White_Noise= np.array([-1.74976547,  0.3426804,   1.1530358,  -0.25243604,  0.98132079], dtype="d")
        #printMatrix(allPoints)
        sim = vicsek.DFDGsim(eigenfunction=vicsek.npEigenVectorsAndValues,s=20, whiteNoise=White_Noise, withConst=True)
        
        #printMatrix(mat)
        #for x in sim:
        #    print(sim[x])

        print(f"vec: {vec}")
        print(f"mat: {mat[:,0]}")
        print("Lap: ")
        printMatrix(laplacian)
        pos = len([x for x in White_Noise if x > 0])
        zeroes = len([x for x in White_Noise if x == 0])
        print(f"zereos: {zeroes}, pos: {pos}, neg: {len(White_Noise)-pos}")
        #printMatrix(laplacian)
        #pprint(f"vec: {vec}")
        #printMatrix(mat)
        #pprint(OGmat)
        
        #print("Doing som calcs. ")
        #j = 1000
        #print(f"eigVal: {vec[j]}")
        #print(f"eigVec: {mat[:, j]}")
#
        #print(vec[j] * mat[:,j])
        #print(laplacian @ mat[:, j])
        #print(np.isclose(vec[j] * mat[:,j], laplacian @ mat[:, j]))

    if whatToRun == 3:
        i=1
        #np.random.seed(100)
        h = 3
        WN = np.random.standard_normal(pow(5, h)*4 + 1)
        #WN = np.random.standard_normal(int((pow(3,h+1)*3)/2))
        #print(f"White Noise: {WN}")
        whiteNoiseList = []
        for h in [1]: 
            WN = np.random.standard_normal(pow(5, h)*4 + 1)
            for withConst in [True]:
                for sVals in [[0.001, 0.5, 1, 20]]:
                    print(f"Vm: {h}")

                    print(f"WN: {WN[0:5]}")
                    theOne = simulation(n = h, s = dh/(2*dw) +0.1, fractalConstruction=vicsekSet.vicsek(n=h).pointsAndNeighbourswhat, fractalType="vicsek")
                    #theOne = simulation(n = h, s = dh/(2*dw) +0.1, fractalConstruction=sg.sierpinski(n=h).pointsAndNeighbours, fractalType="sierpinski")
                    #theOne.drawThePretty(sValues=[sVals],whiteNoise=WN, withConst=withConst)
                    #[0.001, 0.5, 1, 20][20, 50, 75, 100]
                    theOne.drawThePretty(sValues=sVals, whiteNoise=WN, withConst=withConst)
                    #plt.gca().set_position([0, 0, 1, 1])
                    #plt.gcf().set_size_inches(12.8, 4.8)
                    plt.savefig(f"c:/Users/chris/gitProjects/GaussianFieldsThesis/code/sim/Final/{theOne.fractalType}Sim_V{h}_{i}_withConst{withConst}-s{str(sVals).replace(".","_")}.png", format="png", dpi=100)
                    plt.close()
                    whiteNoiseList.append(WN)
                    
    #plt.show()
