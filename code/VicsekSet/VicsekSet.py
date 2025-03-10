import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import timeit
import math
from pprint import pprint

def timer(n: int) -> None:
    """Prints the time for each iteration of the makeTriangle Function

    Args:
    -
        n : int
            To what iteration the loop should run

    Returns:
    -
        None
            Prints the time of each iteration
    """
    for i in range(n):
        start = timeit.default_timer()
        vicsek(n=1).makeVicsek(n = i)
        stop = timeit.default_timer()
        print(f"Iterations:  {i}  Time:  {stop - start}")


class vicsek:
        
    def __init__(self, n):
        """ Initializes the Sierpinski class

        Args:
        -
        n : int
        Defines the number of default reccursions should be done
       
        """
        self.n = int(n) 
    

    def findPointsForLineCollection(self, endPoints: list = None, points: list = None, n: int = None, init_length: int = None) -> list:

        if n == None:
            n = self.n
        if init_length == None:
            init_length = 2
        if points == None:
            points = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        if n <= 0:
            if [[np.array([-init_length, 0]), np.array([init_length, 0])]] in points and [[np.array([0, init_length]), np.array([0,-init_length])]] in points:
                return []
            else:
                return [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        else:
            collection = []

        if endPoints == None:
            endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]


        for x in endPoints:
             collection.extend([[y + x*2 for y in lists] for lists in points])   
        
        newEndPoints = [x + 2*x for x in endPoints]
        
        return collection + self.findPointsForLineCollection(newEndPoints, collection + points, n-1, init_length)

    
    def makeVicsek(self, n: int = None):

        if n == None:
            n = self.n
    
    
        init_length = 2
        init_val = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        init_endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

        if n <= 0:
            line_frac = 1
            midpoints = [(0,0)]
        else:   
            line_frac = n

        vicsekSet = self.findPointsForLineCollection(init_endPoints, init_val, n)
        line_collection = [[tuple(x) for x in lists] for lists in vicsekSet] + [[tuple(x) for x in lists] for lists in init_val]
        midpoints = set([tuple((x[0] + x[1]) // 2) for x in vicsekSet] + [(0,0)])

        # Creates points for scatterplot
        point_x = [x[0] for lists in line_collection for x in lists] + [x[0] for x in midpoints]
        point_y = [y[1] for lists in line_collection for y in lists] + [y[1] for y in midpoints]

        # Create subplots
        fig, ax = plt.subplots(figsize=(10,10))
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1/line_frac, capstyle="butt", joinstyle="round"))
        ax.scatter(point_x, point_y, marker="o", s=5)

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
        y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()
    

    def pointsAndNeighbours(self, n: int = None, init_length: int = None):
        
        if n == None:
            n = self.n
        
        if init_length == None:
            init_length = 3

        vicsekSet = self.findPointsForLineCollection(init_length= init_length)
        allPoints =  {tuple(x) for lists in vicsekSet for x in lists} | set([(0, 0)] + [tuple((x[0] + x[1]) // 2) for x in vicsekSet]) 
        allPointsList = list(allPoints)
        allPointsDict = {}

        for i in range(len(allPoints)):
            neighbours = set()
            for j in [(0, init_length), (0, -init_length), (init_length, 0), (-init_length, 0)]:
                potentialNeighbour = tuple(np.add(allPointsList[i], j))
                if  potentialNeighbour in allPoints:
                    neighbours.add(potentialNeighbour)

            allPointsDict[f"x{i}"] = {"pos" : allPointsList[i], "neighbours" : neighbours, "neighboursSet" : set()}

        return allPointsDict

    def laplacianOperatorMatrix(self):

        allPointsDict = self.pointsAndNeighbours()
        return [[-len(allPointsDict[i]["neighbours"]) if i == j else 1 if allPointsDict[j]["pos"] in allPointsDict[i]["neighbours"] else 0 for i in allPointsDict] for j in allPointsDict]

    def eigenVectorsAndValues(self):
        matrix = self.laplacianOperatorMatrix()
        print(np.diag(matrix))
        return np.linalg.eig(matrix)

    def printLaplacianOperatorMatrix(self) -> None:
        for i in self.laplacianOperatorMatrix():
            print(i)
        return None
    

test = vicsek(n=1)
#test.makeVicsek(n=0)
#pprint(test.pointsAndNeighbours())
w, v = test.eigenVectorsAndValues()

print('E-value:', w)
print('E-vector', v)
#test.makeVicsek()
