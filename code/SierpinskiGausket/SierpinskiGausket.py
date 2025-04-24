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
    To actually draw the triangle one only needs the function "makeTriangle".

    Attributes
    ---------
        n how many recursion of the sierpinski triangle is needed.

    """

    def __init__(self, n):
       """ Initializes the Sierpinski class

       Args
       ----
        n : int
            Defines the number of default recursions should be done
       
       """
       self.n = int(n) 


    def findMidTriangle(self, points: list, n: int) -> list:
        """A recursive function that finds all the upside down triangles of the sierpinski triangle. Meaning it finds all the triangle consisting of the midpoints of the previous triangle

        Args
        ----
            points : list 
                A list of points (often a triangle) for which we find the upside down triangle
            
            n : int
                An integer to measure how far we are in the iteration of the function

        Returns
        -------
            list
                A list of all the upside down triangles until n = 0
        """
        
        # Find the midpoints
        midLeft = self.findMidpoint(points[0], points[2])
        midRight= self.findMidpoint(points[1], points[2])
        midBottom = self.findMidpoint(points[0], points[1])

        # Stops at n-2, because we use integer division and midLeft and MidRight are at (1/4 * 2^n, 2^(n-1)) and (3/4 * 2^n, 2^(n-1))
        if n-2 == 0:
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]]
        else:
            # Gives us the points we just found and run the function again. First on left triangle, then right triangle and lastly the top triangle
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]] + self.findMidTriangle([points[0], midBottom, midLeft], n-1) + self.findMidTriangle([midBottom, points[1], midRight], n-1) + self.findMidTriangle([midLeft, midRight, points[2]], n-1) 


    def findMidpoint(self, startPoint: np.ndarray, endPoint: np.ndarray)-> np.ndarray:
        """ Finds midpoint from two points.
        
        Args
        ----
            StartPoint : np.ndarray
                First point

            EndPoint : np.ndarray
                Second point

        Returns
        -------
            np.ndarray """
        return (startPoint + endPoint) // 2
    

    def trianglePoints(self, n: int = None) -> list:
        """ Find the points in the Sierpinski Gasket from integer. Such you don't have to make your own triangle but can start with a default
        
        Args
        ----
            n : int
                Amount of iterations
        
        Returns
        -------
            list
                List of tuples of coordinates of each point
        """
        
        # Setup if None is given as arguments
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


    def makeTriangle(self, n: int = None, markers: bool = False)-> None:
        """Given an amount of iterations gives a plot for the Sierpinski triangle after n iterations. The plot output will be in the scale of 2**n, to avoid float division and to take advantage of interger division. Further the function "findMidTriangle" finds all the middle triangles, to draw as few lines as possible
        
        Args
        ----
            n : int
                Amount of iterations to run through. If no input is given input from __init__ will be used.
        
        Returns
        -------
            None : 
        """

        # Setup if None is given as arguments
        if n == None:
            n = self.n

        # Create subplots
        fig, ax = plt.subplots()

        sierpinski = self.trianglePoints(n)
        # Add all the lines to the figure
        ax.add_collection(LineCollection(sierpinski, colors = "black", linewidths=1/n))

        if markers == True:
            # Creates points for scatterplot
            point_x = [x[0] for lists in sierpinski for x in lists[:3]]
            point_y = [y[1] for lists in sierpinski for y in lists[:3]]
            
            ax.scatter(point_x, point_y, marker="o", s=5, c='k')

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = pow(2, n+1)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -buffer, pow(2, n+1) + buffer
        y_min, y_max = -buffer, pow(2, n+1) + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()
    

    def pointsWithOrientation(self, points: list = None, n: int = None):
        """ A recursive function that finds all the points to the upside down triangles in the Sierpinski Gasket along with the points "orientation". By orientation we mean if the points are on the "right", "left" or bottom of the triangle. This is so the it can be used later to find each points neighbour.
        
        Args
        ----
            points : list 
                A list of points (often a triangle) for which we find the upside down triangle
            
            n : int
                An integer to measure how far we are in the iteration of the function

        Returns
        -------
            list
                A list of tuples on the following form::
                    
                    [((2, 4), "left"), ((8, 0), "right")... ]
        """

        # Setup if None is given as arguments
        if n == None:
            n = self.n
        if points == None:
            points = [(np.array([0,0]), "left"), (np.array([pow(2, n), 0]), "right"), (np.array([pow(2, n-1),pow(2,n)]), "top")]
        
        # Find midpoints and their orientation
        midLeft = (self.findMidpoint(points[0][0], points[2][0]), "left")
        midRight= (self.findMidpoint(points[1][0], points[2][0]), "right")
        midBottom = (self.findMidpoint(points[0][0], points[1][0]), "bottom")

        # Stops at n-2, because we use integer division and midLeft and MidRight are at (1/4 * 2^n, 2^(n-1)) and (3/4 * 2^n, 2^(n-1))
        if n-2 == 0:
            return [midLeft, midRight, midBottom, (np.array([0,0]), "left"), (np.array([pow(2, self.n), 0]), "right"), (np.array([pow(2, self.n -1),pow(2, self.n)]), "top")]
        else:
            # Gives us the points we just found and run the function again. First on left triangle, then right triangle and lastly the top tirangle
            return [midLeft, midRight, midBottom] + self.pointsWithOrientation([points[0], midBottom, midLeft], n-1) + self.pointsWithOrientation([midBottom, points[1], midRight], n-1)+ self.pointsWithOrientation([midLeft, midRight, points[2]], n-1) 
    

    def pointsAndNeighbours(self, n: int = None): 
        """ A method to find the neighbours of the points of the Sierpinski gasket. The naming convention of the points follow whatever order python chooses in sets.

        Args
        ----
            n : int
                Set the level of recursion
        
        Returns
        -------
            dict
                A dictionary of all points, their coordinate position and their neighbours
        """
        # Setup if None is given as arguments
        if n == None:
            n = self.n

        if n <= 1:
            n = 1
            sierpinskiGasket = {(0,0): "left", (2 ** n, 0): "right", (2**(n-1),2 ** n): "top"}
        else:
            sierpinskiGasket = {(int(x[0][0]), int(x[0][1])): x[1] for x in self.pointsWithOrientation(n=n)}
        
        allPointsDict = {}
        # Nested loops yah :))))))))))
        counter = 0
        for i in sierpinskiGasket:
            # check what is the orientation of a points, so we check the correct neighbours
            if sierpinskiGasket[i] == "left":
                neighboursToCheck = [(1, 2), (2, 0), (1, -2), (-1, -2)]
            elif sierpinskiGasket[i] == "right":
                neighboursToCheck = [(-1, 2), (-2, 0), (-1, -2), (1, -2)]
            elif sierpinskiGasket[i] == "bottom":
                neighboursToCheck = [(-2, 0), (-1, 2), (1, 2), (2, 0)]
            elif sierpinskiGasket[i] == "top":
                neighboursToCheck = [(-1, -2), (1, -2)]
            # Find all potential neighbours
            neighbours = {tuple(np.add(i, j)) for j in neighboursToCheck if tuple(np.add(i, j)) in sierpinskiGasket}
            # Names the points as the appear on allpoints
            allPointsDict[i] = {"neighbours" : neighbours, "name": f"x{counter}"}
            counter +=1

        return allPointsDict