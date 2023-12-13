import numpy as np 
import pandas 
import matplotlib.pyplot as plt
from scipy.optimize import minimize as scipy_minimize
import random
from scipy.optimize import linprog

M=[[j for j in range(3)] for i in range(3)]
print(M)

print(np.dot(M, np.ones(3)))

print(np.linspace(0,1,10)[3])

print(np.random.normal(size=3))

a=np.random.normal(size=10)
id=[1,2]
p=a[id]

print(p, a)