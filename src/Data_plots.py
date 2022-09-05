#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 16:02:57 2022

@author: hxps
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
data = pd.read_csv('test.csv',sep=',')
# con = sqlite3.connect('top-adsorbates-relaxed.db')
# df = pd.read_sql_query("SELECT * FROM tracks", con)
df = pd.DataFrame(data)
df.head()
marker_size = 15 
ratio = df[' Ni_ration']
d = df[' delta_G']
plt.scatter(ratio, d,marker_size, c=d)
plt.ylim((-1,1))
cbar = plt.colorbar()
cbar.set_label('delta_G', labelpad=+1)


plt.show()