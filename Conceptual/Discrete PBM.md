# Attempt at learning generalized doubling rules

## Imports
```python
from aguaclara.play import *
import aguaclara.research.floc_model as floc
import pdb
import pickle as pkl
import matplotlib.animation as anim
import sys
import os
import types
```
## Discrete, same time for collision
```python

pop = np.ones(1000000)
frac = 0.5
collisions = 11
for i in range(1,collisions+1):
    for j in range(0,i):
        # pdb.set_trace()
        sizej = int(2**j)
        popj = np.where(pop==sizej)
        coll = int(np.floor(frac*0.5*np.size(popj)/2)*2)
        if coll > 1:
            pop[popj[0][0:coll]] = 2*sizej
            pop[popj[0][-coll:]] = 0
pop = pop[pop!=0]
bins = np.arange(0,np.log2(np.max(pop)))

plt.hist(pop,bins=bins,align='mid')
bin_w = (max(bins) - min(bins)) / (len(bins) - 1)
plt.xticks(np.arange(min(bins)+bin_w/2, max(bins), bin_w), bins)
plt.xlim(bins[0], bins[-1])
plt.show()
