import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np
import timeit
import math
from pprint import pprint

def timer(n: int) -> None:
    """Prints the time for each iteration of the makeTriangle Function

    Args
    ----
        n : int
            To what iteration the loop should run

    Returns
    -------
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
        """ Initializes the vicsek class

        Args
        ----
            n : int
                Defines the number of default reccursions should be done
       
        """
        self.n = int(n) 
    

    def findPointsForLineCollection(self, endPoints: list = None, points: list = None, n: int = None, init_length: int = None) -> list:
        """A method for finding the lines to draw the vicsek set through the lineCollection method. The method finds the next points by "copying" the whole figure on the endpoints. Since if x is an endpoint and y is a point on the figure then y + 2*x is going to be a new unique point on the next iteration.
        
        Args
        ----
            endPoints : list
                A list of endpoints to "attach" the next iteration to
            points : list
                A list of all the points in the set
            n : int
                recursion depth
            init_length : int
                Lenth of the line pieces
        
        Returns 
        -------
            list
                A list of list with each line that needs to be drawn by lineCollection method

        """

        # Setup if none is given as argument
        if n == None:
            n = self.n
        if init_length == None:
            init_length = 3
        if points == None:
            points = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        
        # This is to check if the recursion given is smaller than zero
        if n <= 0:
            # This checks if we are in the first iteration or if we have just reached recursion depth zero
            if (([[np.array([-init_length, 0]), np.array([init_length, 0])]] and [[np.array([0, init_length]), np.array([0,-init_length])]]) in points):
                return []
            else:
                return [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        else:
            collection = []

        # Setup if none is given as argument
        if endPoints == None:
            endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

        # Finds all the points in the current iteration
        for x in endPoints:
            collection.extend([[y + x*2 for y in lists] for lists in points])

        
        #print(f"collection: {collection}, {len(collection)} \n {points}, end: {endPoints}")
        
        # Find the new endpoints
        newEndPoints = [x + 2*x for x in endPoints]
        
        return collection + self.findPointsForLineCollection(newEndPoints, collection + points, n-1, init_length)


    def findPointsForIDunnoYet(self, endPoints: list = None, points: list = None, n: int = None, init_length: int = None) -> set:
            """A method for finding the lines to draw the vicsek set through the lineCollection method. The method finds the next points by "copying" the whole figure on the endpoints. Since if x is an endpoint and y is a point on the figure then y + 2*x is going to be a new unique point on the next iteration.

            Args
            ----
                endPoints : list
                    A list of endpoints to "attach" the next iteration to
                points : list
                    A list of all the points in the set
                n : int
                    recursion depth
                init_length : int
                    Lenth of the line pieces

            Returns 
            -------
                list
                    A list of list with each line that needs to be drawn by lineCollection method

            """

            # Setup if none is given as argument
            if n == None:
                n = self.n
            if init_length == None:
                init_length = 3
            if points == None:
                points = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
            if endPoints == None:
                endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

                    # This is to check if the recursion given is smaller than zero
            if n <= 0:
                return []
   
            # Finds all the points in the current iteration
            
            

            
            # Finds all the points in the current iteration
            
            #print(f"n: {n}")
            ##print(f"endpoint: {endPoints}")
            ##print(f"y: {set(tuple(np.array(y) + np.array(x)*2) for x in endPoints for lists in points for y in lists)}")
            #pprint(f"points: {points}")
            #pprint(f"huh1: {set([set([tuple([tuple(np.array(y) + np.array(x)*2) for y in lists]) for lists in points]) for x in endPoints])}")
            #pprint(f"test: {set(tuple(y) for x in self.findPointsForLineCollection(n=self.n + 1-n) for y in x)}")
            #collection.update(set([tuple([tuple(np.array(y) + np.array(x)*2)for y in lists]) for lists in points for x in endPoints]))

            #print(f"end: {endPoints}, poin: {points}")
            #for x in endPoints:
            #    collection = [[np.array(y) + 2*np.array(x) for y in lists] for lists in points]
            collection= [[y + x*2 for y in lists] for lists in points for x in endPoints]
            #pprint(f"coll: {collection}")
            #pprint(f"test: {self.findPointsForLineCollection(n=self.n +1 -n)}")
            #pprint(f"testall: {len({tuple(x) for lists in self.findPointsForLineCollection(n=self.n +1 -n) for x in lists} | set([(0, 0)] + [tuple((x[0] + x[1]) // 2) for x in self.findPointsForLineCollection(n=self.n +1 -n)]) )}")

            # Find the new endpoints
            newEndPoints = [x*3 for x in endPoints]

            return collection + self.findPointsForIDunnoYet(endPoints=newEndPoints, points=collection + points, n=n-1, init_length=init_length)
    

    def makeVicsek(self, n: int = None, markers : bool = False):
        """ A mehtod to draw the Vicsek set after n iterations. 

        Args
        ----
            n : int
                recursion depth
            markers : bool
                If there should be markers or not on the plot
        
        Returns
        -------
            None
        
        """

        # Setup if None is given
        if n == None:
            n = self.n
    

        # Set inital values
        init_length = 2
        init_val = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        init_endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

        # Find the points and converts them to tuples so matPlotLib can read them
        vicsekSet = self.findPointsForLineCollection(init_endPoints, init_val, n, init_length)
        line_collection = [[tuple(x) for x in lists] for lists in vicsekSet] + [[tuple(x) for x in lists] for lists in init_val]

        # Setup so the line width is correct
        if n <= 0:
            line_frac = 1
        else:   
            line_frac = n

        # Create subplots
        fig, ax = plt.subplots(figsize=(10,10))
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1/2))

        if markers == True:
            # Creates points for scatterplot
            midpoints = set([tuple((x[0] + x[1]) // 2) for x in vicsekSet] + [(0,0)])
            point_x = [x[0] for lists in line_collection for x in lists] + [x[0] for x in midpoints]
            point_y = [y[1] for lists in line_collection for y in lists] + [y[1] for y in midpoints]
            
            ax.scatter(point_x, point_y, marker="o", s=1, c="k")

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
        y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)

        #plt.grid()
        plt.show()
    

    def pointsAndNeighbourswhat(self, n: int = None, init_length: int = None):
        """A method for getting a dictionary of all the points and neigbors in the Vicsek set. The neighbours are found by calculating all the potential neighbors and checking if they are a point in the Vicsek set. The naming of the points is "random", since the order of the points are decided by the order that they are in the unordered set. This could probably run faster if we had a good naming convention for the points.
        
        Args
        ----
            n : int
                recursion depth
            init_length : int
                Length of the line pieces between any two points
        
        Returns
        -------
            dictionary
                a dictionary of all points who are each a dictionary containing their position and neighbours  
        """

        # Setup if None is given as argument
        if n == None:
            n = self.n        
        if init_length == None:
            init_length = 3

        # Find all the points of the Vicsek set and turn them into a set, such we can do fast lookup.

        if n <= 0:
            allPoints = {(-init_length, 0), (init_length, 0), (0, init_length), (0,-init_length), (0,0)}
        else:
            vicsekSet = self.findPointsForIDunnoYet(n = n, init_length= init_length)
            allPoints =  {tuple(x) for lists in vicsekSet for x in lists} | set([tuple((np.array(x[0]) + np.array(x[1])) // 2) for x in vicsekSet] + [(0, 0)])

        allPointsDict = {} # Initialize dictionary

        # Loop over all points
        for i in allPoints:
            # Find all potential neighbours
            neighbours = {tuple(np.add(i, j)) for j in [(0, init_length), (0, -init_length), (init_length, 0), (-init_length, 0)] if tuple(np.add(i, j)) in allPoints}
            # Creates the dictionary for x_i
            allPointsDict[i] = {"name": f"x{i}", "neighbours" : neighbours}

        return allPointsDict
    

    def laplacianOperatorMatrix(self):
        """Finds the laplacain operators matrix. It is given by L = A - D, where

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

        allPointsDict = self.pointsAndNeighbours()
        return [[-len(allPointsDict[i]["neighbours"]) if i == j else 1 if j in allPointsDict[i]["neighbours"] else 0 for i in allPointsDict] for j in allPointsDict]

    def eigenVectorsAndValues(self):
        matrix = self.laplacianOperatorMatrix()
        return np.linalg.eigh(matrix)

    def printLaplacianOperatorMatrix(self) -> None:
        """ Function to print out the laplacian in a managable way.
        
        Args
        ----
            n : int
                Recursion depth
        """
        for i in self.laplacianOperatorMatrix():
            print(i)
        return None