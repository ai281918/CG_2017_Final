import math
import random
import numpy as np
import matplotlib.tri as tri
from matplotlib import pyplot as plt 

MAX_DISTORTION = 1.5
AREA_LENGTH_RATIO = 0.0

def read_obj(path):
    v = []
    f = []
    obj = open(path)
    lines = obj.readlines()
    for line in lines:
        line = line.split()
        if line[0] == 'v':
            v.append([float(line[1]), float(line[2]), float(line[3])])
        elif line[0] == 'f':
            if len(line) == 4:
                f.append([int(line[1].split('/')[0])-1, int(line[2].split('/')[0])-1, int(line[3].split('/')[0])-1])
            else:
                f.append([int(line[1].split('/')[0])-1, int(line[2].split('/')[0])-1, int(line[4].split('/')[0])-1])
                f.append([int(line[2].split('/')[0])-1, int(line[3].split('/')[0])-1, int(line[4].split('/')[0])-1])
    return np.array(v), np.array(f)

def cross(p1, p2, p3):
    vx1 = p2[0] - p1[0]
    vy1 = p2[1] - p1[1]
    vx2 = p3[0] - p1[0]
    vy2 = p3[1] - p1[1]
    return vx1 * vy2 - vx2 * vy1

def circle_cross(p1, r1, p2, r2):
    # print(p1)
    # print(r1)
    # print(p2)
    # print(r2)
    p3 = np.zeros((2))
    p4 = np.zeros((2))
    # print('abs : ' + str(abs(p1[1] - p2[1])))
    if abs(p1[1] - p2[1]) > 1e-8:
        m = float(p1[0] - p2[0]) / (p2[1] - p1[1])
        k = float(r1*r1 - r2*r2 + p2[0]*p2[0] - p1[0]*p1[0] + p2[1]*p2[1] - p1[1]*p1[1]) / (2 * (p2[1] - p1[1]))
        a = 1 + m*m
        b = 2 * (m*k - m*p2[1] - p2[0])
        c = p2[0]*p2[0] + p2[1]*p2[1] + k*k - 2*k*p2[1] - r2*r2
        if (b*b - 4*a*c) > 0:
            p3[0] = float(-b + math.sqrt(b*b - 4*a*c)) / (2*a)
            p4[0] = float(-b - math.sqrt(b*b - 4*a*c)) / (2*a)
            p3[1] = m*p3[0] + k
            p4[1] = m*p4[0] + k
        elif (b*b - 4*a*c) == 0:
            print('b*b - 4*a*c == 0')
        else:
            print('b*b - 4*a*c < 0')
    else:
        p3[0] = p4[0] = float(-(p1[0]*p1[0] - p2[0]*p2[0] - r1*r1 + r2*r2) / (2*p2[0] - 2*p1[0]))
        a = 1
        b = -2*p1[1]
        c = p3[0]*p3[0] + p1[0]*p1[0] - 2*p1[0]*p3[0] + p1[1]*p1[1] - r1*r1
        p3[1] = float(-b + math.sqrt(b*b - 4*a*c)) / (2*a)
        p4[1] = float(-b - math.sqrt(b*b - 4*a*c)) / (2*a)

    if cross(p1, p2, p3) < 0:
        return p3
    else:
        return p4
    
def plot_triangle(id, f, position_2D, color='blue'):
    p1 = position_2D[f[id][0]]
    p2 = position_2D[f[id][1]]
    p3 = position_2D[f[id][2]]
    # print(p1)
    # print(p2)
    # print(p3)

    rand_data = np.array([[p1[0], p1[1]], [p2[0], p2[1]], [p3[0], p3[1]]])
    triangulation = tri.Triangulation(rand_data[:,0], rand_data[:,1])
    plt.triplot(triangulation, color='blue')

def triangle_area(p1, p2, p3):
    return abs(((p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]))/2)

def length_area_ratio(chart_area, boundary_length, front, id, f, position_2D):
    front = front.split('_')
    new_area = 0
    new_area = triangle_area(position_2D[f[id][0]], position_2D[f[id][1]], position_2D[f[id][2]])
    new_length = 1
    for i in range(3):
        # print(str(f[id][i]) + " " + str(f[id][(i+1)%3]))
        new_length += np.linalg.norm(position_2D[f[id][i]] - position_2D[f[id][(i+1)%3]])
    # print(front)
    new_length -= np.linalg.norm(position_2D[int(front[0])] - position_2D[int(front[1])]) * 2

    # print(str(new_area) + ' ' + str(new_length))
    return (chart_area + new_area) / (boundary_length + new_length), new_area, new_length

def distortiom_metric(id, v, f, position_2D):
    p1 = position_2D[f[id][0]]
    p2 = position_2D[f[id][1]]
    p3 = position_2D[f[id][2]]
    q1 = v[f[id][0]]
    q2 = v[f[id][1]]
    q3 = v[f[id][2]]

    # print('p1 : ' + str(p1))
    # print('p2 : ' + str(p2))
    # print('p3 : ' + str(p3))
    # print('q1 : ' + str(q1))
    # print('q2 : ' + str(q2))
    # print('q3 : ' + str(q3))
    # print((2*triangle_area(p1,p2,p3)))
    # print((2*triangle_area(p1,p2,p3)))
    Ss = (q1*(p2[0]-p3[0])+q2*(p3[0]-p1[0])+q3*(p1[0]-p2[0]))/(2*triangle_area(p1,p2,p3))
    St = (q1*(p3[1]-p2[1])+q2*(p1[1]-p3[1])+q3*(p2[1]-p1[1]))/(2*triangle_area(p1,p2,p3))
    # print('Ss : ' + str(Ss))
    # print('St : ' + str(St))
    a = np.dot(Ss, Ss)
    b = np.dot(Ss, St)
    c = np.dot(St, St)
    # print('a : ' + str(a))
    # print('b : ' + str(b))
    # print('c : ' + str(c))
    Gmax = math.sqrt(((a+c)+math.sqrt((a-c)*(a-c)+4*b*b))/2)
    t = ((a+c)-math.sqrt((a-c)*(a-c)+4*b*b))/2
    if t <= 0:
        return 99999
    Gmin = math.sqrt(t)
    return max(Gmax, 1 / Gmin)

def creat_dictionary(v, f):
    vertex_dic = {}
    edge_dic = {}
    for i, tri in enumerate(f):
        for j in range(3):
            edge = np.array([tri[(j+1)%3], tri[(j+2)%3]])

            if tri[j] not in vertex_dic:
                vertex_dic[tri[j]] = []
            vertex_dic[tri[j]].append(np.concatenate((np.array([i]), edge)))
            # print(vertex_dic[tri[j]])
            
            edge = str(edge[0]) + '_' + str(edge[1])
            edge_dic[edge] = np.array([i, tri[j]])
            # print(edge_dic[edge])
    return vertex_dic, edge_dic
            
def add_front(id, vis, f, fronts, edge_dic):
    for i in range(3):
        edge = str(f[id][(i-1+3)%3]) + '_' + str(f[id][(i-2+3)%3])
        if edge in edge_dic and vis[edge_dic[edge][0]] == 0:
            fronts.add(edge)

def calculate_position(id, v, f, position_2D):
    for i in range(3):
        if f[id][i] not in position_2D:
            # print('Q1 : ' + str(v[f[id][(i+1)%3]]))
            # print('Q2 : ' + str(v[f[id][(i+2)%3]]))
            # print('Q3 : ' + str(v[f[id][i]]))
            # print('R1 : ' + str(np.linalg.norm(v[f[id][i]]-v[f[id][(i+1)%3]])))
            # print('R2 : ' + str(np.linalg.norm(v[f[id][i]]-v[f[id][(i+2)%3]])))
            position_2D[f[id][i]] = circle_cross(position_2D[f[id][(i+1)%3]], np.linalg.norm(v[f[id][i]]-v[f[id][(i+1)%3]]), position_2D[f[id][(i+2)%3]], np.linalg.norm(v[f[id][i]]-v[f[id][(i+2)%3]]))
            # print('P1 : ' + str(position_2D[f[id][(i+1)%3]]))
            # print('P2 : ' + str(position_2D[f[id][(i+2)%3]]))
            # print('P3 : ' + str(position_2D[f[id][i]]))

def insertion():
    return False

def chart_growth(seed_id, v, f, ori, vis, vertex_dic, edge_dic):
    fronts = set([])
    position_2D = {}
    # put the seed triangle
    vis[seed_id] = 1
    add_front(seed_id, vis, f, fronts, edge_dic)
    position_2D[f[seed_id][0]] = np.zeros(2) + ori
    position_2D[f[seed_id][1]] = np.array([0, np.linalg.norm(v[f[seed_id][0]]-v[f[seed_id][1]])]) + ori
    calculate_position(seed_id, v, f, position_2D)
    if seed_id <= 10:
        plot_triangle(seed_id, f, position_2D)
    chart_area = triangle_area(position_2D[f[seed_id][0]], position_2D[f[seed_id][1]], position_2D[f[seed_id][2]])
    boundary_length = 0
    for i in range(3):
        boundary_length += np.linalg.norm(position_2D[f[seed_id][i]] - position_2D[f[seed_id][(i+1)%3]])
    # return
    cnt = 1
    while len(fronts) != 0:
        # print(len(fronts))
        flag = True
        add_list = []
        for front in fronts:
            id = edge_dic[front][0]
            calculate_position(id, v, f, position_2D)
            ratio, new_area, new_length = length_area_ratio(chart_area, boundary_length, front, id, f, position_2D)
            if vis[id] == 0 and distortiom_metric(id, v, f, position_2D) < MAX_DISTORTION and ratio > AREA_LENGTH_RATIO and not insertion():
                # print(chart_area, boundary_length)
                # print(ratio)
                vis[id] = 1
                add_list.append(id)
                if seed_id <= 10:
                    plot_triangle(id, f, position_2D)
                chart_area += new_area
                boundary_length += new_length
                flag = False
        fronts = set([])
        for x in add_list:
            add_front(x, vis, f, fronts, edge_dic)
        cnt += 1
        # if cnt == 1000:
        #     break

def mesh_parameterization(v, f):
    vis = np.zeros(len(f))
    vertex_dic, edge_dic = creat_dictionary(v, f)
    cnt = 0
    for i in range(len(f)):
        if vis[i] == 0:
            print(i)
            chart_growth(i, v, f, np.array([cnt*4, 0]), vis, vertex_dic, edge_dic)
            cnt += 1
            # if cnt == 2:
            #     break
        

    

v, f = read_obj('QQ.obj')
print('Triangle num : ' + str(len(f)))
# print(len(f))
mesh_parameterization(v, f)
# q1 = v[f[0][0]]
# q2 = v[f[0][1]]
# q3 = v[f[0][2]]
# print(q1)
# print(q2)
# print(q3)
# p1 = np.array([0, 0])
# p2 = np.array([0, np.linalg.norm(q1-q2)])
# p3 = circle_cross(p1, np.linalg.norm(q1-q3), p2, np.linalg.norm(q2-q3))
# plot_triangle(p1, p2, p3)
# print(circle_cross(np.array([0,0]), 2.828427, np.array([2,0]), 2))
plt.show()