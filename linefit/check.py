#!/usr/bin/python3

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

df = pd.read_csv('corrected.csv', delimiter=',',header=None,names=['x','y'])

#plt.plot(df['x'], df['y'], label='')
#
#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Title')
#plt.show()
slope, intercept, r_value, p_value, std_err = linregress(df['x'], df['y'])

print ("Slope: " + str(slope))
print ("Intercept: " + str(intercept))

