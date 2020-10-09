#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 15 15:18:53 2020

@author: bolelang
"""

import json
import numpy as np
from pyproj import Proj, transform

# read json data from file

points = []

with open("../frontend/geostreamViz/geostreamvis/gsvApp/public_html/enl_data.json", "r") as read_file:
    el_data = json.load(read_file)
    
    # print(el_data)
    for i in range(len(el_data["data"])):
        coords= el_data["data"][i]["featureOfInterest"]["geometry"]["coordinates"]
        points.append(coords)
        
    #transform form wgs84 to google mercator
    points=np.array(points)
    gm_points = transform(Proj('epsg:4326'), Proj('epsg:3857'), points[:,0], points[:,1])
    print(gm_points)