import networkx as nx
import matplotlib.pyplot as plt
import SierpinskiGausket.SierpinskiGausket as sg
G = nx.Graph()

def what(n):

    if n == 1:    
        sgTri = sg.sierpinski(n=2)
        points = sgTri.trianglePoints()
        i = 0
        for x in points:
            for y in x:
                G.add_node(i, pos=y)
                i += 1

    elif n == 2:
        G.add_node(1, pos=(1,2))
        G.add_node(2, pos=(2,4))
        G.add_node(3, pos=(1,3))
        
    else:
        print("Nothing here buckaroo") 


what(2)
nx.draw(G)
plt.show()