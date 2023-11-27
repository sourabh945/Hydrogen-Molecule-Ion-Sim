from matplotlib import pyplot as plt 
import numpy as np
from mpl_toolkits import mplot3d 
from math import sinh , cosh , sqrt 
import pandas as pd
import os 
from tqdm import tqdm


data = []
e = []
radial = []
n = []
angular = []
d = 2
with open('radial.dat','r') as file : 
    for line in file : 
        data += [line.split()]
        e.append(float(data[-1][0]))
        radial.append(float(data[-1][1]))
    file.close()

with open('angular.dat','r') as file :
    for line in file:
        data += [line.split()]
        n.append(float(data[-1][0]))
        angular.append(float(data[-1][1]))
    file.close()


del data
x = []
y = []
si= []
si_2 = []
data  = []

pbar = tqdm(total=len(e)*len(n),ascii=False,colour='green',ncols=100,desc='Fetching the data :: ',unit="data")

for i in range(0,len(e)):
    for j in range(0,len(n)):
        x.append((d/2)*cosh(e[i])*n[j])
        y.append((d/2)*sinh(e[i])*(sqrt(1-n[j]**2)))
        si.append(radial[i]*angular[j])
        si_2.append(si[i]*si[i])
        data.append([x[-1],y[-1],si[-1],si_2[-1]])
        pbar.update(i*j)

pbar.close()

print(data)

x = np.array(x)
y = np.array(y)
SI = np.array(radial)*np.array(angular)
print(SI)


X,Y = np.meshgrid(x,y)

header = ['x','y','si','si^2']   
data = pd.DataFrame(data,columns=header)

os.system('touch "final_result.csv"')

data.to_csv('final_result.csv',index=False)

print("All data is written in file ")

print("Plotting the figure")

fig = plt.figure(figsize=(15,15))
ax = plt.axes(projection='3d')

surface = ax.plot_surface(X,Y,SI)

ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('z', labelpad=20)

fig.colorbar(surf, shrink=0.5, aspect=8)

plt.show()
