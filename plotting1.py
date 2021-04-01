import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d import axes3d
from matplotlib import cm
import matplotlib.colors as mcolors
import matplotlib.patches as patches
import numpy as np
import mpl_toolkits.mplot3d.art3d as art3d

#ecc = 0
file_name = "3d_results_21_-30.dat"
with open(file_name) as fin:
    fin.readline()
    l = fin.readline()
    pos_start = l.find('eccentricity= ')
    pos_start += len('eccentricity= ')
    pos_end = l.find(' Skew=')
    # print(pos_start, pos_end)
    ecc_str = l[pos_start:pos_end]
    # print(ecc_str)
    ecc = float(ecc_str)
    
    pos_start1 = l.find('Skew= ')
    pos_start1 += len('Skew= ')
    pos_end1 = 101
    # print(pos_start, pos_end)
    skew_str = l[pos_start1:pos_end1]
    # print(ecc_str)
    skew = float(skew_str)
    
 
skew = skew * np.pi/180  


# exit()
X = pd.read_csv(file_name, sep="\t",header = None , skiprows= 8)
X_2d = pd.read_csv("results_c.dat", sep=" ",header = None , skiprows= 8, nrows = 101)

cosphi = np.tile(X_2d[6],101)
sinphi = np.tile(X_2d[7],101)

x = np.zeros_like(X[0])
y = np.zeros_like(X[0])
z = np.zeros_like(X[0])
x_int = np.zeros(100)
y_int1 = np.zeros(100)
z_int = np.zeros(100)
v1n = np.zeros_like(X[0])
v1t = np.zeros_like(X[0])
v2n = np.zeros_like(X[0])
v2t = np.zeros_like(X[0])
sigma1 = np.zeros_like(X[0])
sigma2 = np.zeros_like(X[0])
sigma_int1 = np.zeros(100)
d1l = np.zeros_like(X[0])
d2l = np.zeros_like(X[0])
dl1x = np.zeros_like(X[0])
dl1y = np.zeros_like(X[0])
dl2x = np.zeros_like(X[0])
dl2y = np.zeros_like(X[0])
x = X[0]
y = X[1]
z = X[2]
v1n = X[3]
v1t = X[4]
v2n = X[5]
v2t = X[6]
sigma1 = X[7]
sigma2 = X[8]
d1l = X[11]
d2l = X[12]











#for horizontal lines z = 0 ranges are up [0:101:1] and down [5059:5160:1]
#for vertical lines ranges are up [2534:2635:1] and down [7584:7685:1]
row_start_up = 0
row_end_up = 101
row_start_down = 5059
row_end_down = 5160

# first cooling layer
dl1x = d1l * np.cos(skew) 
dl1y = d1l * np.sin(skew) 
#dl1x = d1l * np.cos(skew) *sinphi
#dl1y = d1l * np.sin(skew) *cosphi
xl1 = x - dl1x
yl1 = y - dl1y
zl1 = z
x1_z0_1 = xl1[row_start_up:row_end_up]
x2_z0_1 = xl1[row_start_down:row_end_down]
y1_z0_1 = yl1[row_start_up:row_end_up]
y2_z0_1 = yl1[row_start_down:row_end_down]
z1_z0_1 = zl1[row_start_up:row_end_up]
z2_z0_1 = zl1[row_start_down:row_end_down]

#second cooling layer
dl2x = d2l * np.cos(skew) 
dl2y = d2l * np.sin(skew)
#dl2x = d2l * np.cos(skew) *sinphi
#dl2y = d2l * np.sin(skew) *cosphi
xl2 = x + dl2x
yl2 = y + dl2y 
zl2 = z
x1_z0_2 = xl2[row_start_up:row_end_up]
x2_z0_2 = xl2[row_start_down:row_end_down]
y1_z0_2 = yl2[row_start_up:row_end_up]
y2_z0_2 = yl2[row_start_down:row_end_down]
z1_z0_2 = zl2[row_start_up:row_end_up]
z2_z0_2 = zl2[row_start_down:row_end_down]


x1_z0 = x[row_start_up:row_end_up]
x2_z0 = x[row_start_down:row_end_down]
y1_z0 = y[row_start_up:row_end_up]
y2_z0 = y[row_start_down:row_end_down]
z1_z0 = z[row_start_up:row_end_up]
z2_z0 = z[row_start_down:row_end_down]

#for widths 
xl_conc = np.concatenate((xl1, xl2))
yl_conc = np.concatenate((yl1, yl2))
zl_conc = np.concatenate((zl1, zl2))
x_z0_l = np.concatenate((x1_z0_1,x2_z0_1,x1_z0_2,x2_z0_2))
y_z0_l = np.concatenate((y1_z0_1,y2_z0_1,y1_z0_2,y2_z0_2))
z_z0_l = np.concatenate((z1_z0_1,z2_z0_1,z1_z0_2,z2_z0_2))

#for other parameters
parameter1_z0  = v1t[row_start_up:row_end_up]
parameter2_z0  = v1t[row_start_down:row_end_down:1]

parameter_z0 = np.concatenate((parameter1_z0,parameter2_z0))
x_z0 = np.concatenate((x1_z0,x2_z0))
y_z0 = np.concatenate((y1_z0,y2_z0))
z_z0 = np.concatenate((z1_z0,z2_z0))

sizes_z0 = np.zeros(len(parameter_z0))
sizes_z0 = 4
sizes = np.zeros(len(X[7]))
sizes = 0.2

fig1 = plt.figure()
axis_lim = 5


ax = fig1.add_subplot(111, projection='3d', aspect = 'auto')
ax.view_init(elev = None, azim = None)

plt.xlabel('x')
plt.ylabel('y')
ax.set_xlim(-axis_lim-1, axis_lim+1)
ax.set_ylim(-axis_lim-1, axis_lim+1)
ax.set_zlim(-axis_lim, axis_lim)

#norm = mcolors.TwoSlopeNorm(vmin = sigma1.min(), vmax= sigma1.max(), vcenter =0.03)
img = ax.scatter(0, 0, 0, c = 'black')
img = ax.scatter(1, 0, 0, c = 'black')
img = ax.scatter(x, y, z, s=sizes, c= d1l, cmap= 'jet')
#img = ax.scatter(x_z0, y_z0, z_z0, s= sizes_z0, c= parameter_z0, cmap= 'jet')








# xy lines for cooling layers
#img = ax.scatter(x_z0_l, y_z0_l, z_z0_l, s= sizes, c = 'red', alpha = 0.3)  
#img = ax.scatter(x_z0, y_z0, z_z0, s= sizes, c = 'black')




plt.colorbar(img)

width = 0.5
height = ((1 - ecc**2) * width**2)**0.5
p = patches.Ellipse((0.5, 0), 2*width, 2*height, fill=False)
ax.add_patch(p)
art3d.pathpatch_2d_to_3d(p, z=0, zdir="z")


#max d1l
x_z0_full = np.concatenate((x1_z0,x2_z0,x1_z0,x2_z0))
y_z0_full = np.concatenate((y1_z0,y2_z0,y1_z0,y2_z0))


l1_for_comparison = ((x_z0_full - x_z0_l)**2 + (y_z0_full - y_z0_l)**2)**0.5
l1max = max(l1_for_comparison)
i_max = l1_for_comparison.argmax()
xmax_values= [x_z0_full[i_max],x_z0_l[i_max]]
ymax_values = [y_z0_full[i_max],y_z0_l[i_max]]

#max d2l
x2l_comp = np.concatenate((x1_z0_2,x2_z0_2))
y2l_comp = np.concatenate((y1_z0_2,y2_z0_2))

l2_for_comparison = ((x_z0 - x2l_comp)**2 + (y_z0 - y2l_comp )**2)**0.5
l2max = max(l2_for_comparison)
j_max = l2_for_comparison.argmax()
xmax2_values = [x_z0[j_max], x2l_comp[j_max]]
ymax2_values = [y_z0[j_max], y2l_comp[j_max]]
print(l2max,j_max)


fig2 = plt.figure()
ax1 = fig2.add_subplot(111, aspect = 'equal')

ax1.set_xlim(-axis_lim-1, axis_lim+1)
ax1.set_ylim(-axis_lim-1, axis_lim+1)

plt.xlabel('x')
plt.ylabel('y')
#img1 = ax1.scatter(x_z0_l, y_z0_l, s = 1, c = 'red')
img1 = ax1.scatter(x_z0, y_z0, s = sizes*2, c = 'black')
img1 = ax1.scatter([0,1],[0,0], c = 'black')
#img1 = ax1.plot([x_z0_full, x_z0_l], [y_z0_full, y_z0_l], lw = 1.5 ,c = 'red', alpha = 0.4)

for i in range(1,101):
    img1 = ax1.plot([x_z0_full[i], x_z0_l[i]], [y_z0_full[i], y_z0_l[i]], lw = 1.5 ,c = 'red', alpha = 0.4)
for i in range(101,202):
    img1 = ax1.plot([x_z0_full[i], x_z0_l[i]], [y_z0_full[i], y_z0_l[i]], lw = 1.5 ,c = 'red', alpha = 0.4)
for i in range(203,303):
    img1 = ax1.plot([x_z0_full[i], x_z0_l[i]], [y_z0_full[i], y_z0_l[i]], lw = 1.5 ,c = 'red', alpha = 0.4)
for i in range(303,404):
    img1 = ax1.plot([x_z0_full[i], x_z0_l[i]], [y_z0_full[i], y_z0_l[i]], lw = 1.5 ,c = 'red', alpha = 0.4)

img1 = ax1.plot(xmax2_values, ymax2_values, lw = 1.5 ,c = 'green', alpha = 1)

img1 = ax1.plot(xmax_values, ymax_values, lw = 1.5 ,c = 'green', alpha = 1)


plt.axhline(y=0.0,c = 'black', lw = sizes*2)
#plt.axvline(x=x[0])

plt.show()
