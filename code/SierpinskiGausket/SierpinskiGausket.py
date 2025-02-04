import matplotlib.pyplot as plt
import numpy as np

class sierpinski:
    def __init__(self, n):
       self.n = int(n) 


    def findMidTriangle(self, points: list, n: int):
        
        midLeft = self.findMidPoint()

        if n == 0:
            self.draw(points)
        else:

        

    def draw(self, points: list):
        x = np.array([x[0] for x in points]+[points[0][0]])
        y = np.array([y[1] for y in points]+[points[0][1]])
        plt.plot(x, y)
    
    def findMidpoint(self, startPoint: np.ndarray, endPoint):
        return (startPoint + endPoint) // 2

    def makeTriangle(n: int):
        init_val = [np.array([0,0]), np.array([2 ** n, 0]), np.array([2**(n-1),2 ** n])]
        self.draw(init_val)



        plt.grid()
        plt.show()
        


test = sierpinski(n=10)
n = 1
init_val = [np.array([0,0]), np.array([2 ** n, 0]), np.array([2**(n-1),2 ** n])]

test.draw(init_val)