import numpy as np

def find_adjacent(coord, resolution):
    xs = []
    if coord[0] > 0:
        xs.append(coord - np.array([1, 0]))
    if coord[0] < resolution[0] - 1:
        xs.append(coord + np.array([1, 0]))
    if coord[1] > 0:
        xs.append(coord - np.array([0, 1]))
    if coord[1] < resolution[1] - 1:
        xs.append(coord + np.array([0, 1]))
    return xs

def solve_eikonal(coord, resolution, u, f):
    if coord[0] == 0:
        xmin = u[coord[0] + 1, coord[1]]
    elif coord[0] == resolution[0] - 1:
        xmin = u[coord[0] - 1, coord[1]]
    else:
        xmin = np.min([u[coord[0] + 1, coord[1]], u[coord[0] - 1, coord[1]]])
    if coord[1] == 0:
        ymin = u[coord[0], coord[1] + 1]
    elif coord[1] == resolution[0] - 1:
        ymin = u[coord[0], coord[1] - 1]
    else:
        ymin = np.min([u[coord[0], coord[1] + 1], u[coord[0], coord[1] - 1]])
    ans = np.min([xmin, ymin]) + 1 / f[coord[0], coord[1]]
    if not ans > np.max([xmin, ymin]):
        return ans
    else:
        ans = (xmin + ymin + np.sqrt(- (xmin - ymin) ** 2 + 2 / f[coord[0], coord[1]] ** 2)) / 2
        return ans

def put_in_list(inlist, Q, coord):
    for item in coord:
        inlist[item[0], item[1]] = 1
        Q.append(item)

## initialize
source = [0, 0]
L = [16, 16]
epsilon = 10 ** -2
resolution = np.array([16, 16])
inlist = np.zeros(resolution)
dx = L / resolution
speed = np.ones(resolution)
u = np.ones(resolution) * np.Inf
u[source[0], source[1]] = 0
Q = []
put_in_list(inlist, Q, find_adjacent(source, resolution))
## update
while len(Q) > 0:
    temp = Q.pop(0)
    p = u[temp[0], temp[1]]
    q = solve_eikonal(temp, resolution, u, speed)
    u[temp[0], temp[1]] = q
    if np.abs(p - q) < epsilon:
        neighbors = find_adjacent(temp, resolution)
        for item in neighbors:
            if inlist[item[0], item[1]] == 0:
                p = u[item[0], item[1]]
                q = solve_eikonal(item, resolution, u, speed)
                if p > q:
                    u[item[0], item[1]] = q
                    put_in_list(inlist, Q, [item])
        inlist[temp[0], temp[1]] = 0
    else:
        put_in_list(inlist, Q, [temp])
print(u)
