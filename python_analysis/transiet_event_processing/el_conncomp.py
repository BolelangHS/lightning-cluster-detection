#!/usr/bin/env python3

import scipy as sp
import numpy as np
import math
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
import datetime

# to read json data only
import json
from pyproj import Proj, transform

# for alphashapes only
import sys
from descartes import PolygonPatch
import alphashape
from shapely.geometry import Polygon, LineString, Point

# for plotting maps
import geopandas as gpd

""" Here we use a fixed point to repesent an event. An event is not neccesarily a point so this class
needs to be generic enough to handle any count of event 
"""

'''class Point:
  def __init__(self, x, y):
    self.x = x
    self.y = y
    self.id = -1'''
    
def read_points(file_name, points, gm_points, start_date, delta):
    
    end_date = start_date + datetime.timedelta(minutes=delta)
    with open(file_name, "r") as read_file:
        el_data = json.load(read_file)
        
        # print(el_data)
        for i in range(len(el_data["data"])):
            data_time = datetime.datetime.fromisoformat(el_data["data"][i]["phenomenonTime"])
            
            if start_date < data_time and data_time < end_date:
                coords= el_data["data"][i]["featureOfInterest"]["geometry"]["coordinates"]
                points.append(coords)
            
        if len(points) > 0:
            #transform form wgs84 to google mercator
            points=np.array(points)
            x_points, y_points = transform(Proj('epsg:4326'), Proj('epsg:3857'), points[:,0], points[:,1], always_xy=True)
            for i in range(len(x_points)):
                if y_points[i] < -2800000:
                    gm_points.append([x_points[i], y_points[i]])
            #print(x_points[i])
            # print('points .........')
            # print(points)
            # print('gm_points .........')
            # print(gm_points)
        

# Create a delaunay
def generate_delaunay(points):

    points = np.array(points)
    if len(points) >= 3:
        triangulation = Delaunay(points)
        #world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
        #world = world.to_crs("EPSG:3857") 
        #south_africa = world[world.name == "South Africa"]
        sa_provinces = gpd.read_file("/home/bsibolla/Data/sa_provinces/South_Africa_files/sa_provinces_wgs84.shp")
        sa_provinces = sa_provinces.to_crs("EPSG:3857")
        #sa_center_north = sa_provinces[sa_provinces.ID == "Free State"]
        sa_center_north = sa_provinces[(sa_provinces.ID != "Western Cape")&(sa_provinces.ID != "Eastern Cape")]
        sa_center_north.plot(alpha=0.1, edgecolor='brown')
        # ax = south_africa.boundary.plot()
        plt.triplot(points[:,0], points[:,1], triangulation.simplices, color='grey', linewidth=0.5)
        plt.plot(points[:,0], points[:,1], 'o', markersize=4, color='red', alpha=0.3)
        plt.ticklabel_format(style='plain')
        plt.xticks(rotation=45)
        plt.show()
        #print(triangulation.simplices)
        return triangulation

def distance(pointa, pointb):
    dis = math.sqrt(pow(pointb[0] - pointa[0], 2) + pow(pointb[1] - pointa[1], 2))
    return dis
    
# Lets create the graph for connected components    
def create_graph(triangulation, size, points, threshold):
    
    matrix = csr_matrix((size, size), dtype=np.int8).toarray()
    
    # iterate over triangles in the delaunay - to add edges to the graph
    for i in range(len(triangulation.simplices)):
        triangle = triangulation.simplices[i]
        # Get the length of each edge between nodes
        dis_a = distance(points[triangle[0]], points[triangle[1]]);
        dis_b = distance(points[triangle[1]], points[triangle[2]]);
        dis_c = distance(points[triangle[2]], points[triangle[0]]);
        # to print out distances
        # print('dis a...............')
        # print(dis_a)
        # print('dis_b...............')
        # print(dis_b)
        # print('dis_c ..............')
        # print(dis_c)   
        # only adding edges that meet filtering criteria - enforcing filters 
        if dis_a < threshold:
            matrix[triangle[0]][triangle[1]] = 1
        if dis_b < threshold:
            matrix[triangle[1]][triangle[2]] = 1
        if dis_c < threshold:
            matrix[triangle[2]][triangle[0]] = 1
            
    graph = csr_matrix(matrix)
    # print('matrix ...............')
    # print(matrix)   
    return graph  

# get components and plot
def print_components(points, labels, component_size):
    
    # bunch points in to components
    components = {}
    for c in range(component_size):
        components[c] = []
    
    for i in range(len(points)):
        components[labels[i]].append(points[i])
    
    # plot the components
    sa_provinces = gpd.read_file("/home/bsibolla/Data/sa_provinces/South_Africa_files/sa_provinces_wgs84.shp")
    sa_provinces = sa_provinces.to_crs("EPSG:3857")
    sa_center_north = sa_provinces[(sa_provinces.ID != "Western Cape")&(sa_provinces.ID != "Eastern Cape")]
        
    # world = gpd.read_file(gpd.datasets.get_path('naturalearth_lowres'))
    # world = world.to_crs("EPSG:3857") 
    # south_africa = world[world.name == "South Africa"]
    #ax = south_africa.boundary.plot()
    #fig, ax = plt.subplots()
    #set_aspect('equal')
    sa_center_north.plot(alpha=0.1, edgecolor='brown')
    #south_africa.boundary.plot()
    for c in range(component_size):
        comp_points = np.array(components[c])
        plt.ticklabel_format(style='plain')
        plt.xticks(rotation=45)
        plt.plot(comp_points[:,0], comp_points[:,1], 'o', markersize=5, alpha=0.3)
        
    plt.show()
    # print('components ...............')
    # print(components)
    return components
    
def create_alphashapes(components, alpha, time, result):
       
    for i in range(len(components)):

        fig, ax = plt.subplots()
        cluster = {}
        cluster["cluster_label"] = i
        cluster["cluster_points"] = components[i]
        cluster["cluster_time"] = time.isoformat()
        
        if len(components[i]) < 3:
            if len(components[i]) == 1:
                polygon = LineString([(components[i][0][0], components[i][0][1]), (components[i][0][0], components[i][0][1])]).buffer(1000, 16, 30)
            elif len(components[i]) == 2:
                polygon = LineString([(components[i][0][0], components[i][0][1]), (components[i][1][0], components[i][1][1])]).buffer(1000, 16, 30)
            
            cluster["cluster_polygon"] = []
        else:
            
            alpha_shape = alphashape.alphashape(np.array(components[i]))
            
            coords = []
            for p in alpha_shape.exterior.coords:
                #coords.append([p[0], p[1]])
                coords.append((p[0], p[1]))
                
            polygon = Polygon(coords).buffer(300, 16, 30)
            
                
        cluster["cluster_polygon"] = list(polygon.exterior.coords)
        # print('component [i]..............')
        # print(components[i])
        # print('alphashape..............')
        # print(alpha_shape)
        ax.ticklabel_format(style='plain')
        ax.tick_params('x', labelrotation=45)
        ax.scatter(*zip(*components[i]))
        ax.add_patch(PolygonPatch(polygon, alpha=0.2))
            
        result["result"].append(cluster)
     
    # print(result)
    plt.show()
    
 
# Global Result object    


daysdata = []
delta = 5 # in minutes
numdays = 1
start_time = datetime.datetime.fromisoformat("2017-11-28T16:30:00+00:00")
end_time = start_time + datetime.timedelta(days=numdays)
curr_time = start_time

daycount = 0
while daycount < numdays:
    times = []
    curr_time = start_time + datetime.timedelta(days=daycount)
    #end_time = curr_time + datetime.timedelta(days=1)
    while curr_time < end_time:
        times.append(curr_time)
        curr_time += datetime.timedelta(minutes=delta)
        
    daysdata.append(times)
    daycount+=1

print(daysdata)

result = {}
result["source"] = ''
result["type"] = ''
result["method"] = ''
result["result"] = []
    
for day in daysdata:
    for ctime in day:
        print('Searching Data for Date ' + ctime.isoformat())
        
    
        # print('Reading for Date ' + ctime.isoformat())
        # Create an empty list of points     
        points = []    
        gm_points = []
        threshold = 4000.00
        
        read_points("../data/enl_data.json", points, gm_points, ctime, delta)
        
        if len(points) > 0:
            print('Data Found for Date ' + ctime.isoformat())

            
            if len(points) > 0 and len(points) < 3:
                dis = 0
                labels = []
                n_components = 0
                
                # if its one point then there is one component
                if len(points) == 1:
                    labels.append(0)
                    n_components = 1
                # if its two point there could be one or 2 components depending  on distance 
                elif len(points) == 2:    
                     dis  = math.sqrt(pow(gm_points[1][0] - gm_points[0][0], 2) + pow(gm_points[1][1] - gm_points[0][1], 2))
                     #print(gm_points)
                     labels.append(0)
                     labels.append(1)
                     n_components = 2
                
                     if dis < threshold:
                        labels[1] = 0
                        n_components = 1
                    
                components = print_components(gm_points, labels, n_components)
                create_alphashapes(components, 5.0, ctime, result)
            elif len(points) >= 3:
                print('Points found for time period of ' + str(delta) + ' minute(s) : ' + str(len(points)))
                # generate triangulation
                tri = generate_delaunay(gm_points)
                
                #generate graph
                sz = len(gm_points)
                graph = create_graph(tri, sz, gm_points, threshold)
                
                # run connected components - to get the clusters
                n_components, labels = connected_components(csgraph=graph, directed=False, return_labels=True)
                # print('n_components ...............')
                # print(n_components)
                # print('labels ...............')
                # print(labels)
                
                components = print_components(gm_points, labels, n_components)
                create_alphashapes(components, 5.0, ctime, result)
    
print('Writing JSON output for ' +  day[0].isoformat())
filename = '../data/results/enl_concom_result.json'
            
with open(filename, 'w') as file:
    file.write(json.dumps(result)) 

    