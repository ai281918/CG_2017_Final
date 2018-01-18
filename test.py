import math
import random
import numpy as np
import matplotlib.tri as tri
from matplotlib import pyplot as plt 
from functools import cmp_to_key

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
    
def plot_triangle(id, f, position_2D, ori, color='blue'):
    p1 = position_2D[f[id][0]]
    p2 = position_2D[f[id][1]]
    p3 = position_2D[f[id][2]]
    flag1 = True
    flag2 = True
    for i in range(3):
        a = abs(position_2D[f[id][i]][0] - position_2D[f[id][(i+1)%3]][0])
        b = abs(position_2D[f[id][i]][1] - position_2D[f[id][(i+1)%3]][1])
        if a > 1e-3:
            flag1 = False
        if b > 1e-3:
            flag2 = False
        if a < 1e-3 and b < 1e-3:
            return
        
    if flag1 or flag2:
        return
    # print(p1)
    # print(p2)
    # print(p3)

    rand_data = np.array([[p1[0]+ori[0], p1[1]+ori[1]], [p2[0]+ori[0], p2[1]+ori[1]], [p3[0]+ori[0], p3[1]+ori[1]]])
    # print(rand_data)
    triangulation = tri.Triangulation(rand_data[:,0], rand_data[:,1])
    plt.axis('equal')
    plt.triplot(triangulation, color='blue')


def polt_texture(texture, f, position_2D):
    l = []
    for i,tex in enumerate(texture):
        Mx = -1
        mx = 1
        My = -1
        my = 1
        for t in tex:
            for j in f[t]:
                if position_2D[i][j][0] > Mx:
                    Mx = position_2D[i][j][0]
                if position_2D[i][j][0] < mx:
                    mx = position_2D[i][j][0]
                if position_2D[i][j][1] > My:
                    My = position_2D[i][j][1]
                if position_2D[i][j][1] < my:
                    my = position_2D[i][j][1]
        l.append([i, Mx,-mx,My,-my])
    l.sort(key=cmp_to_key(lambda a,b:(b[3]+b[4])-(a[3]+a[4])))

    cur_x = 0
    cur_y = 0
    next_y = 0
    next_x = 0
    for i in range(len(l)):
        # print(l[i])
        if i % int(math.sqrt(len(l))) == 0:
            cur_y += next_y + l[i][4]
            next_y = l[i][3] + 0.1
            cur_x = 0
            next_x = 0
        cur_x += next_x + l[i][2]
        next_x = l[i][1] + 0.1
        for id in texture[l[i][0]]:
            plot_triangle(id, f, position_2D[l[i][0]], [cur_x, cur_y])
            # print(distortiom_metric(id, v, f, position_2D[l[i][0]]))

                
    # texture.sort(key=cmp_to_key(lambda a,b:len(b) - len(a)))
    # print(texture)


def triangle_area(p1, p2, p3):
    return max(abs(((p2[0]-p1[0])*(p3[1]-p1[1])-(p3[0]-p1[0])*(p2[1]-p1[1]))/2), 1e-8)

def length_area_ratio(chart_area, boundary_length, front, id, f, position_2D):
    front = front.split('_')
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

def calculate_position(vertex_id, vertex_dic, v, f, fronts, position_2D, vis):
    tmp_list = []
    distortion_sum = 0
    for x in vertex_dic[vertex_id]:
        edge = str(x[1]) + '_' + str(x[2])
        if edge in fronts:
            position_2D[vertex_id] = circle_cross(position_2D[x[1]], np.linalg.norm(v[vertex_id]-v[x[1]]), position_2D[x[2]], np.linalg.norm(v[vertex_id]-v[x[2]]))
            d = distortiom_metric(x[0], v, f, position_2D)
            distortion_sum += d
            tmp_list.append([position_2D[vertex_id], d])
    
    position_2D[vertex_id] = np.zeros(2)
    for x in tmp_list:
        position_2D[vertex_id] += x[0] * x[1] / distortion_sum
    if len(tmp_list) > 0:
        position_2D[vertex_id] = tmp_list[0][0]

def check_distortion(vertex_id, vertex_dic, v, f, fronts, position_2D, vis):
    for x in vertex_dic[vertex_id]:
        edge = str(x[1]) + '_' + str(x[2])
        if vis[x[0]] == 0 and edge in fronts:
            if distortiom_metric(x[0], v, f, position_2D) > MAX_DISTORTION:
                return False
    return True

def insertion():
    return False

def add_vertex_to_chart(vertex_id, add_list, vertex_dic, fronts, texture):
    pass


def chart_growth(seed_id, v, f, vis, vertex_dic, edge_dic, texture):
    fronts = set([])
    position_2D = {}
    # put the seed triangle
    vis[seed_id] = 1
    add_front(seed_id, vis, f, fronts, edge_dic)
    position_2D[f[seed_id][0]] = np.zeros(2)
    position_2D[f[seed_id][1]] = np.array([0, np.linalg.norm(v[f[seed_id][0]]-v[f[seed_id][1]])])
    position_2D[f[seed_id][2]] = circle_cross(position_2D[f[seed_id][0]], np.linalg.norm(v[f[seed_id][2]]-v[f[seed_id][0]]), position_2D[f[seed_id][1]], np.linalg.norm(v[f[seed_id][2]]-v[f[seed_id][1]]))
    texture[-1].append(seed_id)
    chart_area = triangle_area(position_2D[f[seed_id][0]], position_2D[f[seed_id][1]], position_2D[f[seed_id][2]])
    boundary_length = 0
    for i in range(3):
        boundary_length += np.linalg.norm(position_2D[f[seed_id][i]] - position_2D[f[seed_id][(i+1)%3]])
    # add vertex
    while len(fronts) != 0:
        # print(len(fronts))
        add_list = []
        for front in fronts:
            tri_id = edge_dic[front][0]
            vertex_id = edge_dic[front][1]
            pop_flag = vertex_id not in position_2D
            if pop_flag:
                calculate_position(vertex_id, vertex_dic, v, f, fronts, position_2D, vis)
            ratio, new_area, new_length = length_area_ratio(chart_area, boundary_length, front, tri_id, f, position_2D)
            if vis[tri_id] == 0 and check_distortion(vertex_id, vertex_dic, v, f, fronts, position_2D, vis) and ratio > AREA_LENGTH_RATIO and not insertion():
                for x in vertex_dic[vertex_id]:
                    edge = str(x[1]) + '_' + str(x[2])
                    if vis[x[0]] == 0 and edge in fronts:
                        # print(x[0])
                        vis[x[0]] = 1
                        add_list.append(x[0])
                        texture[-1].append(x[0])
                chart_area += new_area
                boundary_length += new_length
            elif pop_flag:
                # pass
                position_2D.pop(vertex_id, None)
        fronts = set([])
        for x in add_list:
            add_front(x, vis, f, fronts, edge_dic)
    
    return position_2D

def mesh_parameterization(v, f):
    texture = []
    position_2D = []
    vis = np.zeros(len(f))
    vertex_dic, edge_dic = creat_dictionary(v, f)
    cnt = 0
    for i in range(len(f)):
        if vis[i] == 0:
            print(i)
            texture.append([])
            position_2D.append(chart_growth(i, v, f, vis, vertex_dic, edge_dic, texture))
            cnt += 1
            # if cnt == 2:
            #     break
    polt_texture(texture, f, position_2D)
        

    

v, f = read_obj('sphere.obj')
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