#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
from scipy.spatial import cKDTree
from metpy.interpolate.geometry import dist_2
from metpy.interpolate.points import barnes_point
from metpy.interpolate.tools import average_spacing, calc_kappa
import pandas as pd
import math

df = pd.read_csv()
keyname="T2M"
t2m=df[keyname]

# Creating the Hourly Index Array

xp=np.arange(0,len(t2m))
zeros = np.zeros(len(t2m))
xp2=list(zip(xp,zeros))

# Changing Floating-Point Temperatures

zp=t2m

### First iteration at all points        
sim_gridx=xp
sim_gridy = np.zeros(len(sim_gridx))
grid_points = np.column_stack((sim_gridx,sim_gridy))

#Restrict to a Radius of Observations
#Radius of influence 
radius = 30
obs_tree = cKDTree(xp2)
indices = obs_tree.query_ball_point(grid_points, r=radius)
 
#Radius of influence
#Fix for Points Within Reach Indices of Influence Radius
x=[]
y=[]
for item in indices:
    u,v = obs_tree.data[item].T
    x.append(u)
    y.append(v)
#Distance squared for each of the points
barnes_dist0 = []
for index,item in enumerate(indices):
    u=dist_2(sim_gridx[index], sim_gridy[index], x[index], y[index])
    barnes_dist0.append(u)

#Observations for each point

barnes_obs0_wn = []
for item in indices:
    barnes_obs0_wn.append(zp[item]) 
#Filtrar los valores válidos
for index, (dist, values) in enumerate(zip(barnes_dist0,barnes_obs0_wn)):
    distsAndValues = list(zip(dist,values))
    distsAndValues =  [(d,v) for (d,v) in distsAndValues if not(pd.isna(v))]
    res=list(zip(*distsAndValues))
    ndist=res[0]
    nobs=res[1]
    barnes_dist0[index]=np.array(list(ndist))
    barnes_obs0_wn[index]=np.array(list(nobs))
    barnes_obs0=barnes_obs0_wn

#Calculate the kappa coefficient
kappa = calc_kappa(average_spacing(xp2))
#Calculating Barnes Approximations for All Points         
#First iteration 
barnes_val0=[]
for dist, obs in zip(barnes_dist0,barnes_obs0):
    barnes_val0.append(barnes_point(dist, obs, kappa,gamma=1))
#Second iteration on the points without observation      
#Searching for points without observation 
sim_gridx1=sim_gridx
sim_gridx1=np.array(sim_gridx1)       
sim_gridy1 = np.zeros(len(sim_gridx1))
grid_points1 = np.column_stack((sim_gridx1,sim_gridy1))
indices1 = obs_tree.query_ball_point(grid_points1, r=radius)
x1=[]
y1=[]
for item in indices1:
    u,v = obs_tree.data[item].T
    x1.append(u)
    y1.append(v)

# Distances for each point
barnes_dist1 = []
for index,item in enumerate(indices1):
    u=dist_2(sim_gridx1[index], sim_gridy1[index], x1[index], y1[index])
    barnes_dist1.append(u)

barnes_obs1_wn = []
for item in indices1:
    barnes_obs1_wn.append(zp[item]-np.array(barnes_val0)[item])
    
# Filter out values without nan
    for index, (dist, values) in enumerate(zip(barnes_dist1,barnes_obs1_wn)):
    distsAndValues = list(zip(dist,values))
    distsAndValues =  [(d,v) for (d,v) in distsAndValues if not(pd.isna(v))]
    res=list(zip(*distsAndValues))
    ndist=res[0]
    nobs=res[1]
    barnes_dist1[index]=np.array(list(ndist))
    barnes_obs1_wn[index]=np.array(list(nobs))
    barnes_obs1=barnes_obs1_wn

#Barn values with modified weight
barnes_val1_temp=[]
for dist, obs in zip(barnes_dist1,barnes_obs1):
    barnes_val1_temp.append(barnes_point(dist, obs, kappa,gamma=0.2))
    
#Calculate Barn Values for the Second Iteration at Unobserved Points

barnes_val1=[x + y for x, y in zip(barnes_val1_temp, barnes_val0)]


    

