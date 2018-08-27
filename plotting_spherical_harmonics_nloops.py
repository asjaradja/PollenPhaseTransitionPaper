from __future__ import division
import scipy as sci
import scipy.special as sp
import matplotlib
from matplotlib import cm, colors
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import cm, colors
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation
import numpy as np
from scipy.special import sph_harm
import sys

#importing file
cms = sys.argv[1]
parameters = sys.argv[2]

Re_list = []
Im_list = []
parameters_list = []


for line in open(cms):
        list = line.split()
        try:
                Re = float(list[0])
                Im = float(list[1])
        except TypeError:
                print("Type Error!")
                continue
        Re_list.append(Re)
        Im_list.append(Im)

for line in open(parameters):
        list = line.split()
        try:
                value = str(list[0])
                parameter = float(list[1])

        except TypeError:
                print("Type Error!")
                continue
        parameters_list.append(parameter)

        if len(parameters_list)>=5:
                break


a = open(parameters)
lines = a.readlines()
if lines:
        last_line=lines[-1]
H_value = last_line.split()[1]

parameters_list.append(H_value)


phi = np.linspace(0, np.pi, 150)
theta = np.linspace(0, 2*np.pi, 150)
theta_2 = np.linspace(np.pi, 3*np.pi, 150)

phi, theta = np.meshgrid(phi, theta)

l = (int)(len(Re_list)/3)-1

R = (Re_list[0]*sph_harm(0,l-1,theta,phi)).real + (Re_list[l]*sph_harm(0,l,theta,phi)).real + (Re_list[2*l+1]*sph_harm(0,l+1,theta,phi)).real
for n in range(1,l):
        R += ((complex(Re_list[n],Im_list[n])*sph_harm(n,l-1,theta,phi)).real) + ((complex(Re_list[n],Im_list[n]).conjugate()*pow(-1,n)*sph_harm(-n,l-1,theta,phi)).real)
for n in range(1,l+1):
        R += ((complex(Re_list[n+l],Im_list[n+l])*sph_harm(n,l,theta,phi)).real) + ((complex(Re_list[n+l],Im_list[n+l]).conjugate()*pow(-1,n)*sph_harm(-n,l,theta,phi)).real)
for n in range(1,l+2):
	R += ((complex(Re_list[n+2*l+1],Im_list[n+2*l+1])*sph_harm(n,l+1,theta,phi)).real) + ((complex(Re_list[n+2*l+1],Im_list[n+2*l+1]).conjugate()*pow(-1,n)*sph_harm(-n,l+1,theta,phi)).real)

# The Cartesian coordinates of the unit sphere
s = 0.02
X = (s*R-1) * np.sin(phi) * np.cos(theta)
Y = (s*R-1) * np.sin(phi) * np.sin(theta)
Z = (s*R-1) * np.cos(phi)
#second side of sphere
x_2 = (s*R-1) * np.sin(phi) * np.cos(theta_2)
y_2 = (s*R-1) * np.sin(phi) * np.sin(theta_2)

#Plotting Spherical Harmonics
norm = colors.Normalize()
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,2,1, projection='3d')
#ax.plot_surface(X, Y, Z,  rstride=1, cstride=1, color='y',edgecolor='none')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm(R)))
ax.view_init(0,0)
ax.set_axis_off()
ax = fig.add_subplot(1,2,2, projection='3d')
ax.plot_surface(X, Y, Z, rstride=1, cstride=1, facecolors=cm.jet(norm(R)))
ax.view_init(0,180)

ax.set_title('Parameters are: el_not=%s, tau_bar=%s, lambda_bar=%s, E_thresh=%s, num_ls=%s\nMinimum Hamiltonian value: %s'%(parameters_list[0], parameters_list[1], parameters_list[2], parameters_list[3], parameters_list[4], parameters_list[5]), fontsize=8)
ax.set_axis_off()
plt.savefig("el_not_%s_tau_bar_%s_lambda_bar_%s_E_thresh_%s_num_ls_%s.png"%(parameters_list[0], parameters_list[1], parameters_list[2], parameters_list[3], parameters_list[4]))
#plt.show()

#plot mayavi interactive plot
#from mayavi import mlab
#mlab.figure(1, bgcolor=(0,0,0), fgcolor=(0,0,0), size=(400,400))
#mlab.clf()
#mlab.mesh(X,Y,Z,scalars=norm(R),colormap='YlOrBr')
#mlab.view(90.70,6.2, (-1.3,-2.9,0.25))
#mlab.show()
