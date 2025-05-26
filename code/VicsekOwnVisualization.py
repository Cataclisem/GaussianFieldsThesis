import matplotlib.pyplot as plt
import numpy as np

n = 4

def genFunc(point, i):
    #print(f"i: {i}, point: {point}")
    return 1/3*(np.array(point)-np.array(i))+np.array(i)

init_length=pow(3,n)
start_points = [np.array([-init_length, 0]), np.array([init_length, 0]), np.array([0, init_length]), np.array([0,-init_length]), np.array([0,0])]


collection = []
beforePoints = start_points.copy()

if n <= 1:
    collection = start_points.copy()
else:
    for counter in range(1,n):
        for start_point in start_points:
            collection = collection + [genFunc(x, start_point) for x in beforePoints]

        beforePoints += collection

x_points, y_points = zip(*[(x[0], x[1]) for x in collection])

reset = []
for start_point in start_points:
    #theNext.extend([genFunc(x, start_point) for x in beforePoints])
    #print(f"before: {beforePoints}, len: {len(beforePoints)}")
    some = [genFunc(x, start_point) for x in beforePoints]
    #print(f"some: {some}, len: {len(some)}")
    reset.extend(some)
    #print(f"before2: {beforePoints}, len: {len(beforePoints)}")
    #a = input("Stopping: ")

print(f"lenres: {len(reset)}")

collection.extend(reset)
theNext = []
for start_point, col in zip(start_points, ["r", "g", "b", "y", "m"]):
    #theNext.extend([genFunc(x, start_point) for x in beforePoints])
    x_pointsNext, y_pointsNext = zip(*[(x[0], x[1]) for x in [genFunc(x, start_point) for x in collection]])
    plt.scatter(x_pointsNext, y_pointsNext, marker='o', color=col, s=2)

#print(f"x: {x_points}, \n y: {y_points}")

print(f"lenres2: {len(reset)}")

x_pointsReset, y_pointsReset = zip(*[(x[0], x[1]) for x in reset])
plt.scatter(x_pointsReset, y_pointsReset, marker='o', color="k", s=4)

plt.scatter(x_points, y_points, marker='o', color="tab:brown", s=8)


plt.show()