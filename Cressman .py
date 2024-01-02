#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# Calculate the distances between each point and all other points
distances = np.zeros((len(dfmi), len(dfmi)))
for i in range(len(dfmi)):
    for j in range(len(dfmi)):
        distances[i,j] = np.sqrt((dfmi['Dia'][i] - dfmi['Dia'][j])**2 + (dfmi['Hora'][i] - dfmi['Hora'][j])**2)

# Radius of influence in units of length

radius = 10  # 
weights = np.zeros((len(dfmi), len(dfmi)))
for i in range(len(dfmi)):
    for j in range(len(dfmi)):
        if i != j and distances[i,j] <= radius:
            weights[i,j] = 1 - (distances[i,j] / radius)

# Normalize the weights for each point
for i in range(len(dfmi)):
    weight_sum = np.sum(weights[i,:])
    if weight_sum > 0:
        weights[i,:] /= weight_sum

# Calculates the interpolated values for each missing point in the column
for i in range(len(dfmi)):
    if np.isnan(dfmi['T2M'][i]):
        z_interpolated = 0
        weight_sum = 0
        for j in range(len(dfmi)):
            if i != j and not np.isnan(dfmi['T2M'][j]):
                z_interpolated += weights[j,i] * dfmi['T2M'][j]
                weight_sum += weights[j,i]
        if weight_sum > 0:
            dfmi['T2M'][i] = z_interpolated / weight_sum


