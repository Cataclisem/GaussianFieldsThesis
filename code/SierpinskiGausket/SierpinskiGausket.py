import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from pprint import pprint
import timeit

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
        sierpinski(n=1).makeTriangle(n = i)
        stop = timeit.default_timer()
        print(f"Iterations:  {i}  Time:  {stop - start}")


class sierpinski:
    """
    A class to contain all the functions to draw the sierpinski triangle. 
    To actually draw the trinangle one only needs the function "makeTriangle".

    Attributes:
    -
        n : int
            How many recursion of the sierpinski triangle is needed.

    """

    def __init__(self, n):
       """ Initializes the Sierpinski class

       Args:
       -
        n : int
            Defines the number of default reccursions should be done
       
       """
       self.n = int(n) 


    def findMidTriangle(self, points: list, n: int) -> list:
        """A recursive function that finds all the "inner" triangles of the sierpinski triangle.

        Meaning it finds all the triangle consisting of the midpoints of the previous triangle

        Args:
        -
            points : list 
                A list of points (often a triangle) for which we find the "inner" triangle
            
            n : int
                An integer to measure how far we are in the iteration of the function

        Returns:
        -   
            list
                A list of all the inner triangles until n = 0
        """
        
        # Find the midpoints
        midLeft = self.findMidpoint(points[0], points[2])
        midRight= self.findMidpoint(points[1], points[2])
        midBottom = self.findMidpoint(points[0], points[1])

        # Stops at n-2, because we use interger division and midLeft and MidRight are at (1/4 * 2^n, 2^(n-1)) and (3/4 * 2^n, 2^(n-1))
        if n-2 == 0:
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]]
        else:
            # Gives us the points we just found and run the function again. First on left triangle, then right triangle and lastly the top tirangle
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]] + self.findMidTriangle([points[0], midBottom, midLeft], n-1) + self.findMidTriangle([midBottom, points[1], midRight], n-1) + self.findMidTriangle([midLeft, midRight, points[2]], n-1) 


    def pointsWithOrientation(self, points: list = None, n: int = None, i: int = None):
        if n == None:
            n = self.n
        if points == None:
            points = [(np.array([0,0]), "left"), (np.array([2 ** n, 0]), "right"), (np.array([2**(n-1),2 ** n]), "top")]
        midLeft = (self.findMidpoint(points[0][0], points[2][0]), "left")
        midRight= (self.findMidpoint(points[1][0], points[2][0]), "right")
        midBottom = (self.findMidpoint(points[0][0], points[1][0]), "bottom")
        if n-2 == 0:
            return [midLeft, midRight, midBottom]
        else:
            # Gives us the points we just found and run the function again. First on left triangle, then right triangle and lastly the top tirangle
            return [midLeft, midRight, midBottom] + self.pointsWithOrientation([points[0], midBottom, midLeft], n-1) + self.pointsWithOrientation([midBottom, points[1], midRight], n-1)+ self.pointsWithOrientation([midLeft, midRight, points[2]], n-1) 
    
    def pointsAndNeighbours(self, n: int = None): 
        if n == None:
            n = self.n
        
        sierpinskiGasket = {(tuple(x[0]), x[1]) for x in self.pointsWithOrientation(n = n)} | {(tuple(x[0]), x[1]) for x in [(np.array([0,0]), "left"), (np.array([2 ** n, 0]), "right"), (np.array([2**(n-1),2 ** n]), "top")]}
        allPoints, orientation = zip(*sierpinskiGasket)
        allPointsList = list(allPoints)
        allPointsDict = {}

        for i in range(len(allPoints)):
            neighbours = set()
            if orientation[i] == "left":
                neighborsToCheck = [(1, 2), (2, 0), (1, -2), (-1, -2)]
            elif orientation[i] == "right":
                neighborsToCheck = [(-1, 2), (-2, 0), (-1, -2), (1, -2)]
            elif orientation[i] == "bottom":
                neighborsToCheck = [(-2, 0), (-1, 2), (1, 2), (2, 0)]
            elif orientation[i] == "top":
                neighborsToCheck = [(-1, -2), (1, -2)]

            for j in neighborsToCheck:
                potentialNeighbour = tuple(np.add(allPointsList[i], j))
                if  potentialNeighbour in allPoints:
                    neighbours.add(potentialNeighbour)

            allPointsDict[f"x{i}"] = {"pos" : allPointsList[i], "neighbours" : neighbours, "neighboursSet" : set()}

        return allPointsDict
    
    def laplacianOperatorMatrix(self):

        allPointsDict = self.pointsAndNeighbours()
        return [[-len(allPointsDict[i]["neighbours"]) if i == j else 1 if allPointsDict[j]["pos"] in allPointsDict[i]["neighbours"] else 0 for i in allPointsDict] for j in allPointsDict]


    def printLaplacianOperatorMatrix(self) -> None:
        for i in self.laplacianOperatorMatrix():
            print(i)
        return None

    def eigenVectorsAndValues(self):
        matrix = self.laplacianOperatorMatrix()
        print(np.diag(matrix))
        return np.linalg.eig(matrix)


    def findMidpoint(self, startPoint: np.ndarray, endPoint: np.ndarray)-> np.ndarray:
        """ Finds midpoint from two points.
        
        Args:
        -
            StartPoint : np.ndarray
                First point

            EndPoint : np.ndarray
                Second point

        Returns:
        --
            np.ndarray """
        return (startPoint + endPoint) // 2
    
    def trianglePoints(self, n: int = None) -> list:
        
        # If no n choosen just choose the one define from the class
        if n == None:
            n = self.n
        
        if n <= 0:
            n = 1
            return [[tuple(x) for x in [np.array([0,0]), np.array([2, 0]), np.array([1,2]), np.array([0,0])]]]
        else:
            n = n + 1 
             # Sets up initial Triangle
            init_val = [np.array([0,0]), np.array([2 ** n, 0]), np.array([2**(n-1),2 ** n]), np.array([0,0])]
            return self.findMidTriangle(init_val, n) + [[tuple(x) for x in init_val]]

    def makeTriangle(self, n: int = None)-> None:
        """Given an amount of iterations gives a plot for the Sierpinski triangle after n iterations.
            The plot output will be in the scale of 2**n, to avoid float division and to take advantage of interger division

            Futher the function "findMidTriangle" finds all the middle triangles, to draw as few lines as possible
        
        Args:
        -
            n : int (default = None)
                Amount of iterations to run through. If no input is given input from __init__ will be used.
        
        Returns:
        -
            None
        """

        if n == None:
            n = self.n
        

        line_collection = self.trianglePoints(n)

        # Create subplots
        fig, ax = plt.subplots()
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1/n))

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = pow(2, n)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -buffer, pow(2, n) + buffer
        y_min, y_max = -buffer, pow(2, n) + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()


test = sierpinski(n=2)

test.printLaplacianOperatorMatrix()