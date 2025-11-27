import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import sum, label
M = 2000
Llist = [50, 200, 450]
a = 1.2
plist = [0.5, 0.52, 0.54, 0.56, 0.58]
allarea = np.array([])

for r in range(1, len(Llist)+1):
    plt.figure(figsize=(8,8))
    L = Llist[r-1]
    for d in range(1, len(plist)+1):
        allarea = np.array([])
        p = plist[d-1]
        for i in range(M):
            z = np.random.rand(L,L)
            m = z<p
            lw, num = label(m)
            labelList = np.arange(lw.max() + 1)
            area = sum(m,lw,labelList)
            allarea = np.append(allarea,area)
        n,sbins = np.histogram(allarea,bins=int(max(allarea)))
        s = 0.5*(sbins[1:]+sbins[:-1])
        logmax = np.ceil(np.log(max(s))/np.log(a))
        logbins = a**np.arange(0,logmax)
        nl,nlbins = np.histogram(allarea,bins=logbins)
        ds=np.diff(logbins)
        sl=0.5*(logbins[1:]+logbins[:-1])
        nsl=nl/(M*L**2*ds)
        #plt.subplot(3,1,r)
  
        plt.loglog(sl,nsl,'.', label=f'p = {p}')

    plt.xlabel('$s$')
    plt.ylabel('$n(s,p;L)$')
    plt.title(f'$L = {L}$, $a = 1.2$')   
    plt.legend()   
plt.show()
