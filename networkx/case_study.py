#!/usr/bin/env python
# coding: utf-8

# In[1]:


import networkx as nx
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')


# In[2]:


Alpha = nx.read_edgelist("datasets/bitcoinalpha.csv", delimiter = ',', create_using=nx.DiGraph, data=(('weight', int),('sec', float),))


# In[3]:


Ce = nx.closure(Alpha)
Ce_weighted = nx.closure(Alpha, weight = 'weight')


# In[4]:


degree_ce = dict()
for k, v in Ce.items():
    #if v > 0:    #remove 0, in order to use log
    degree_ce[k] = [Alpha.degree(k), v[0]]
print(len(degree_ce))
        
# 1 and -1 to assign color in scatter plot
degree_ce_weight = dict()   
for k, v in Ce_weighted.items():
    if v[0] >= 0:     
        degree_ce_weight[k] = [Alpha.degree(k, 'weight'), v[0], 1]  # delete/add 'weight' to compare with degree
    else:
        degree_ce_weight[k] = [Alpha.degree(k, 'weight'), v[0], -1]  # delete/add 'weight' to compare with degree
print(len(degree_ce_weight))

degree_ce_weight_negative = dict()   
for k, v in Ce_weighted.items():
    if v[0] < 0:
        degree_ce_weight_negative[k] = [Alpha.degree(k, 'weight'), v[0]]   # delete/add 'weight' to compare with degree
print(len(degree_ce_weight_negative))

degree_ce_weight_high = dict()   
for k, v in Ce_weighted.items():
    if v[0] > 0.25:
        degree_ce_weight_high[k] = [Alpha.degree(k, 'weight'), v[0]]   # delete/add 'weight' to compare with degree
print(len(degree_ce_weight_high))


# In[5]:


print(nx.average_normalized_patterns_app(Alpha,degree_ce_weight_negative.keys()))
print(nx.average_normalized_patterns_app(Alpha,degree_ce_weight_high.keys()))


# In[6]:


degree_ce_list = []
for key, value in degree_ce.items():
    temp = [key, value[0], value[1]]
    degree_ce_list.append(temp)
degree_ce_list = sorted(degree_ce_list, key=lambda t:t[2])

degree_ce_weight_list = []
for k, v in degree_ce_weight.items():
    temp = [k, v[0], v[1], v[2]]
    degree_ce_weight_list.append(temp)
degree_ce_weight_list = sorted(degree_ce_weight_list, key=lambda t:t[2], reverse=True)

degree_ce_negative_list = []
for key, value in degree_ce_weight_negative.items():
    temp = [key, value[0], value[1]]
    degree_ce_negative_list.append(temp)
degree_ce_negative_list = sorted(degree_ce_negative_list, key=lambda t:t[2])

#degree_ce_weight_high = degree_ce_weight_list[0:15]
degree_ce_high_list = []
for key, value in degree_ce_weight_high.items():
    temp = [key, value[0], value[1]]
    degree_ce_high_list.append(temp)
degree_ce_high_list = sorted(degree_ce_high_list, key=lambda t:t[2])


# In[7]:


df_degree_ce = pd.DataFrame(degree_ce_list, columns=['node-id','degree', 'clo'])
df_degree_ce_weighted = pd.DataFrame(degree_ce_weight_list, columns=['node-id','strength', 'clo', 'sign'])
df_degree_ce_negative = pd.DataFrame(degree_ce_negative_list, columns=['node-id','strength', 'clo'])
df_degree_ce_high = pd.DataFrame(degree_ce_high_list, columns=['node-id','strength', 'clo'])


# In[8]:

#Fig 1
onepic, axes = plt.subplots(nrows=2,ncols=1,figsize=(8,12.5))
plt.subplots_adjust(hspace=0.25)  
fig = df_degree_ce.plot(kind="scatter", x = 'degree', y = 'clo',ax=axes[0], logx=True, s=1, c='Black')
fig.set_xlabel("degree", fontsize=14)
fig.set_ylabel("clo. coef.", fontsize=14)
fig.set_title("BTC-Alpha (w/o weight)", fontsize=16)
fig.set_xlim(0.9, 1200)
fig.set_ylim(-0.001, 0.6)
fig.tick_params(axis = 'x', labelsize = 14)
fig.tick_params(axis = 'y', labelsize = 14)
fig.annotate("ρ = 0.714", xy = (16, 0.52), xytext=(16, 0.52), fontsize = 16)

fig_weight = df_degree_ce_weighted.plot(kind = 'scatter', x = 'strength', y = 'clo', ax=axes[1], s=1, c='sign', colormap='coolwarm', colorbar=False)
fig_weight.set_xlabel("strength", fontsize=14)
fig_weight.set_ylabel("clo. coef.", fontsize=14)
fig_weight.tick_params(axis = 'x', labelsize = 14)
fig_weight.tick_params(axis = 'y', labelsize = 14)
fig_weight.set_title("BTC-Alpha (w/ weight)", fontsize=16)
fig_weight.annotate("ρ = 0.265", xy = (20, 0.50), xytext=(20, 0.50), fontsize = 16)
plt.savefig('alpha.pdf')


#Fig 2
pic, axes = plt.subplots(nrows=1,ncols=2,figsize=(13,4))
fig_negative = df_degree_ce_negative.plot(kind = 'scatter', ax=axes[0], x = 'strength', y = 'clo', s=5, c='blue')
fig_negative.set_xlabel("strength", fontsize=14)
fig_negative.set_ylabel("clo. coef.", fontsize=14)
fig_negative.set_title("BTC-Alpha (clo. coef. < 0)", fontsize=16)
fig_negative.tick_params(axis = 'x', labelsize = 14)
fig_negative.tick_params(axis = 'y', labelsize = 14)

fig_high = df_degree_ce_high.plot(kind='scatter', ax=axes[1], x = 'strength', y = 'clo', s=10, c='red')
fig_high.set_xlabel("strength", fontsize=14)
fig_high.set_ylabel("", fontsize=14)
fig_high.set_title("BTC-Alpha (clo. coef. > 0.3)", fontsize=16)
fig_high.tick_params(axis = 'x', labelsize = 14)
fig_high.tick_params(axis = 'y', labelsize = 14)
plt.savefig('alpha_zoom.pdf')


# Three observations
# 1. highly trusted group are only formed with small circle
# 2. nodes with big strength, is better then average, but not as much as when weight is ignored
# 3. very small number of negative triangle, and not too low.





