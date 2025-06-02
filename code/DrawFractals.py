import SierpinskiGausket.SierpinskiGausket as sg
import VicsekSet.VicsekSet as vicsekSet
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from matplotlib.collections import LineCollection
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import numpy as np
from pprint import pprint



vicsek = vicsekSet.vicsek(n = 0)

def makeVicsek(n: int = vicsek.n, markers: bool = False):    

        # Set inital values
        init_length = 2
        init_val = [[np.array([-init_length, 0]), np.array([init_length, 0])], [np.array([0, init_length]), np.array([0,-init_length])]]
        init_endPoints = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length])]

        # Find the points and converts them to tuples so matPlotLib can read them
        vicsekSet = vicsek.findPointsForLineCollection(init_endPoints, init_val, n, init_length)
        line_collection = [[tuple(x) for x in lists] for lists in vicsekSet] + [[tuple(x) for x in lists] for lists in init_val]

        # Setup so the line width is correct
        if n <= 0:
            line_frac = 1
        else:   
            line_frac = n

        # Create subplots
        fig, ax = plt.subplots(figsize=(10,10))
        ax.add_collection(LineCollection(line_collection, colors = "black", linewidths=1))

        if markers == True:
            # Creates points for scatterplot
            midpoints = set([tuple((x[0] + x[1]) // 2) for x in vicsekSet] + [(0,0)])
            point_x = [x[0] for lists in line_collection for x in lists] + [x[0] for x in midpoints]
            point_y = [y[1] for lists in line_collection for y in lists] + [y[1] for y in midpoints]
            
            ax.scatter(point_x, point_y, marker="o", s=0.2, c="k")

        # Computes limites of graph (how far x and y axis should stretch out) based on n
        buffer = (pow(3, n) * init_length)/10 # A buffer to make the final graph not look as cramped
        x_min, x_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer
        y_min, y_max = -pow(3, n) * init_length - buffer, pow(3, n) * init_length + buffer

        ax.set_xlim(x_min, x_max)
        ax.set_ylim(y_min, y_max)


        plt.axis('off')
        plt.gca().set_position([0, 0, 1, 1])
        plt.savefig(f"c:/Users/chris/gitProjects/GaussianFieldsThesis/code/vicsek{n}Markers{str(markers)}.svg", format="svg")
        plt.close(fig)
        #plt.show()

makeVicsek(n = 3, markers=True)