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
    def __init__(self, n):
       self.n = int(n) 


    def findMidTriangle(self, points: list, n: int) -> list:
        
        midLeft = self.findMidpoint(points[0], points[2])
        midRight= self.findMidpoint(points[1], points[2])
        midBottom = self.findMidpoint(points[0], points[1])

        # Stops at n-2, because we use interger division and midLeft and MidRight are at 1/4 * 2^n
        if n-2 == 0:
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]]
        else:
            return [[tuple(midLeft), tuple(midRight), tuple(midBottom), tuple(midLeft)]] + self.findMidTriangle([points[0], midBottom, midLeft], n-1) + self.findMidTriangle([midBottom, points[1], midRight], n-1) + self.findMidTriangle([midLeft, midRight, points[2]], n-1) 
            
        
    
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

    def makeTriangle(self, n: int = None)-> None:
        """Given an amount of iterations gives a plot for the Sierpinski triangle after n iterations.
            The plot output will be in the scale of 2**n, to avoid float division and to take advantage of interger division
        
        Args:
        -
            n : int
                Amount of iterations to run through 
        
        Returns:
        -
            None
        """

        # If no n choosen just choose the one define from the class
        if n == None:
            n = self.n
        n = n + 1 
        
        # Sets up initial Triangle
        init_val = [np.array([0,0]), np.array([2 ** n, 0]), np.array([2**(n-1),2 ** n]), np.array([0,0])]
        
        # If initial n <= 0, then we just draw the triangle
        if n <= 1:
            line_collection = [[tuple(x) for x in init_val]]
        else:
            line_collection = self.findMidTriangle(init_val, n) + [[tuple(x) for x in init_val]]

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
        


test = sierpinski(n=7)
test.makeTriangle()