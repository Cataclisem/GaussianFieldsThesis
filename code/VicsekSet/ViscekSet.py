import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


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
        
        if n == 0:
            return collection + points

        newEndPoints = [x + 2*x for x in endPoints]
        return collection + self.findPoints(newEndPoints, collection + points, n-1)

    
    def makeViscek(self, n: int = None):

        if n == None:
            n = self.n


        init_length = 2
        init_val = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        init_endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

        test = self.findPoints(init_endPoints, init_val, n)

        line_collection = [[tuple(x) for x in lists] for lists in test]
        #print(test)

        # Create subplots
        fig, ax = plt.subplots()
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1/n))

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
        y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()

        
        return None

test = viscek(n=4)
test.makeViscek()