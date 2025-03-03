import numpy as np
import sys

sys.path.append("../PyCompareC/")
from CompareC import *

np.random.seed(42)

rej_perc = []
for _ in range(2000):
    print(_)
    nn = 100
    lifetimes = np.random.exponential(scale=np.exp(1), size=nn)
    censtimes = np.random.exponential(scale=10, size=nn)
    x1 = np.random.normal(loc=lifetimes)
    x2 = np.random.normal(loc=lifetimes)
    ztimes = np.minimum(lifetimes, censtimes)
    status = (censtimes > lifetimes).astype(int)
    
    rej_perc.append(compareC(ztimes, status, x1, x2)['pval'])

print(np.mean(np.array(rej_perc) < 0.05))  # Expected ~0.0465
print(compareC(ztimes, status, x1, x2))
