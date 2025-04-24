import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import numpy as np


class hataTree:

    def __init__(self, n: int):
        self.n = n
    

    def recursion(self, points: list = None, n: int = None, alpha: np.array = np.array([-0.4, 0.65])):
        if n == None:
            n = self.n
        if points == None:
            points = [np.array([0,0]), np.array([0,1])]
        
        if n <= 0:
            return []

        print(f"n: {n}")

        len_alpha = np.linalg.norm(alpha)

        collectionF1 = [np.array([alpha[0]*x[0] + alpha[1]*x[1], alpha[1]*x[0] - alpha[0]*x[1]]) for x in points]
        collectionF2 = [pow(len_alpha, 2)+(1-pow(len_alpha, 2))*x*np.array([1,-1]) for x in points]

        return collectionF1 + collectionF2 + self.recursion(points= collectionF1+collectionF2+points, n=n-1, alpha=alpha)
    

    def recursionWithLines(self, points: list = None, n: int = None, alpha: np.array = np.array([-0.4, 0.65])):
        if n == None:
            n = self.n
        if points == None:
            points = [(np.array([0,0]), np.array([0,1])), (np.array([0,1]), np.array([0,0])), (alpha, np.array([0,0]))]
        
        if n <= 0:
            return []

        print(f"n: {n}")

        len_alpha = np.linalg.norm(alpha)

        collectionF1 = [[np.array([alpha[0]*x[0][0] + alpha[1]*x[0][1], alpha[1]*x[0][0] - alpha[0]*x[0][1]]), x[0]] for x in points]
        collectionF2 = [[pow(len_alpha, 2)+(1-pow(len_alpha, 2))*x[0]*np.array([1,-1]), x[0]] for x in points]

        return collectionF1 + collectionF2 + self.recursionWithLines(points= collectionF1+collectionF2+points, n=n-1, alpha=alpha)
    

    def plot(self, points):
        
        pointsx, pointsy = zip(*[(x[0], x[1]) for x in points])
        plt.scatter(pointsx, pointsy, s=1)
        plt.show()

    
    def plotLines(self, points):
        
        # Create subplots
        fig, ax = plt.subplots()

        # Add all the lines to the figure
        ax.add_collection(LineCollection(points, colors = "black", linewidths=1/self.n))

        plt.show()

test = hataTree(n=1)

print(test.recursionWithLines())
test.plotLines(test.recursionWithLines())
#test.plot(test.recursion())