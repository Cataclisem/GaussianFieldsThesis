import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import timeit
import math

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
        viscek(n=1).makeViscek(n = i)
        stop = timeit.default_timer()
        print(f"Iterations:  {i}  Time:  {stop - start}")


class viscek:
        
    def __init__(self, n):
        """ Initializes the Sierpinski class

        Args:
        -
        n : int
        Defines the number of default reccursions should be done
       
        """
        self.n = int(n) 
    

    def findPoints(self, endPoints: list, points: list, n: int = None) -> list:

        if n == None:
            n = self.n
    
        collection = []

        for x in endPoints:
             collection.extend([[y + x*2 for y in lists] for lists in points])
        
        if n-1 == 0:
            return collection

        newEndPoints = [x + 2*x for x in endPoints]
        return collection + self.findPoints(newEndPoints, collection + points, n-1)

    
    def makeViscek(self, n: int = None):

        if n == None:
            n = self.n


        init_length = 2
        init_val = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]

        if n <= 0:
            line_frac = 1
            line_collection = init_val
            midpoints = [(0,0)]
        else:   
            line_frac = n
            init_endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]
            vicsekSet = self.findPoints(init_endPoints, init_val, n)
            line_collection = [[tuple(x) for x in lists] for lists in vicsekSet] + [[tuple(x) for x in lists] for lists in init_val]
            midpoints = set([tuple((x[0] + x[1]) // 2) for x in vicsekSet] + [(0,0)])
        #print(test)
        
        point_x = [x[0] for lists in line_collection for x in lists] + [x[0] for x in midpoints]
        point_y = [y[1] for lists in line_collection for y in lists] + [y[1] for y in midpoints]

        # Create subplots
        fig, ax = plt.subplots(figsize=(10,10))
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1/line_frac, capstyle="butt", joinstyle="round"))
        ax.scatter(point_x, point_y, marker="o", s=20)

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
        y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()

viscek(n=1).makeViscek(n=0)
#timer(10)