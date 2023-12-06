from matplotlib import pyplot as plt 
import csv
import numpy as np
from math import pow
import os 
from mpl_toolkits import mplot3d

''' this program we use to operate on the data file we generate through 
the fortran program here we are plot the graph with help of this program 
here we make all the plots except some of the are made using the gnuplot 
'''
path = "./processed_csv_files"

e = []
n = []
angular = []
radial = []

data = []

os.system(f"touch {path}/radial.csv")
os.system(f"touch {path}/angular.csv")
os.system(f"{path}/Hydrogen_Radial.csv")


with open('radial.dat','r') as file :
    with open(f'{path}/radial.csv','a') as csv_file:
        writer = csv.writer(csv_file)
        for line in file : 
            data += [line.split()]
            e.append(float(data[-1][0]))
            radial.append(float(data[-1][1]))
            writer.writerow([e[-1],radial[-1]])
        csv_file.close()
    file.close()

with open('angular.dat','r') as file :
    with open(f'{path}/angular.csv','a') as csv_file:
        writer = csv.writer(csv_file)
        for line in file:
            data += [line.split()]
            n.append(float(data[-1][0]))
            angular.append(float(data[-1][1]))
            writer.writerow([n[-1],angular[-1]])
        csv_file.close()
    file.close()

Radial = []
E = []
step = 20/200000

with open(f'{path}/Hydrogen_Radial.csv',"a") as file :
    writer = csv.writer(file)
    for i in range(len(e)-1,-1,-1):
        Radial.append(radial[i])
        writer.writerow([radial[i]])
    for i in range(0,len(e)):
        Radial.append(radial[i])
        writer.writerow([radial[i]])
    
    file.close()

result = []
separation = 200

with open(f'{path}/HMI_Radial.csv','a') as file :
    writer = csv.writer(file)
    for i in range (0,len(Radial)+separation):
        if i > separation and i < len(Radial):
            result.append(Radial[i]+Radial[i-separation])
        elif i>=len(Radial):
            result.append(Radial[i-separation])
        else:
            result.append(Radial[i])
        writer.writerow([result[i]])

resultab = []

with open(f'{path}/HMI_Radial_antibonding.csv','a') as file :
    writer = csv.writer(file)
    for i in range (0,len(Radial)+separation):
        if i > separation and i < len(Radial):
            resultab.append(-Radial[i]+Radial[i-separation])
        elif i>=len(Radial):
            resultab.append(Radial[i-separation])
        else:
            resultab.append(-Radial[i])
        writer.writerow([resultab[i]])
        
Result = np.array(result)

Resultab = np.array(resultab)


plt.plot(Radial)
plt.title('Hydrogen Radial Function from HMI Radial function')
plt.ylabel('R(e)')
plt.xlabel('Arbitrary axis')
plt.savefig(f"{path}/Hydrogen Radial Function")
plt.close()
plt.plot(result)
plt.title(f'HMI Radial function at nuclear axis separation {separation}')
plt.ylabel('R(e)')
plt.xlabel('Arbitrary axis')
plt.savefig(f"{path}/HMI Radial Function{separation}")
plt.close()
plt.plot(resultab)
plt.title(f'HMI Radial function at nuclear axis anti bonding separation {separation}')
plt.ylabel('R(e)')
plt.xlabel('Arbitrary axis')
plt.savefig(f"{path}/HMI Radial Function ab {separation}")
plt.close()

step = 10/len(Result)
step_y = 10/len(Radial)
x = np.arange(0,10,step_y)
y = np.arange(0,10,step)
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(12,10))
ax = plt.axes(projection='3d')

Result = np.array(Result)
print('multiply start')
Z = []
with open(f'{path}/HMI_Radial.csv','a') as file :
    writer = csv.writer(file)
    for i in range(0,len(Result)):
        mid = []
        for j in range(0,len(Radial)):
            mid.append((Result[i]*Radial[j]))
        writer.writerow(mid)
        Z.append(mid)

Z = np.array(Z)
print('complete')

surf = ax.plot_surface(X,Y,Z,cmap=plt.cm.cividis)

plt.title('HMI bonding state')

ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('R(e)', labelpad=20)

fig.colorbar(surf, shrink=0.5, aspect=8)

plt.savefig(f"HMI bonding state {separation*10}.png")

plt.show(block=True)

step = 10/len(Resultab)
step_y = 10/len(Radial)
x = np.arange(0,10,step_y)
y = np.arange(0,10,step)
X,Y = np.meshgrid(x,y)

fig = plt.figure(figsize=(12,10))
ax = plt.axes(projection='3d')

Resultab= np.array(Resultab)
print('multiply start')
Z = []
with open(f'{path}/HMI_Radial.csv','a') as file :
    writer = csv.writer(file)
    for i in range(0,len(Resultab)):
        mid = []
        for j in range(0,len(Radial)):
            mid.append((Resultab[i]*Radial[j]))
        writer.writerow(mid)
        Z.append(mid)

Z = np.array(Z)
print('complete')

surf = ax.plot_surface(X,Y,Z,cmap=plt.cm.cividis)

plt.title('HMI anti bonding state')

ax.set_xlabel('x', labelpad=20)
ax.set_ylabel('y', labelpad=20)
ax.set_zlabel('R(e)', labelpad=20)

fig.colorbar(surf, shrink=0.5, aspect=8)

plt.savefig(f"HMI anti bonding state {separation*10}.png")

plt.show(block=True)