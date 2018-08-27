import matplotlib as mpl
#mpl.use('Agg')
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

	if len(parameters_list)>=4:
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

# The Cartesian coordinates of the unit sphere
x = np.sin(phi) * np.cos(theta)
y = np.sin(phi) * np.sin(theta)
z = np.cos(phi)

x_2 = np.sin(phi) * np.cos(theta_2)
y_2 = np.sin(phi) * np.sin(theta_2)

#m and l named correctly
l = len(Re_list)-1
# Calculate the spherical harmonic Y(l,m) and normalize to [0,1]
fcolors = (Re_list[0]*sph_harm(0,l,theta,phi)).real
for m in range(1,l+1):
	fcolors += ((complex(Re_list[m],Im_list[m])*sph_harm(m,l,theta,phi)).real) + ((complex(Re_list[m],Im_list[m]).conjugate()*pow(-1,m)*sph_harm(-m,l,theta,phi)).real) 	

fcolors_2 = (Re_list[0]*sph_harm(0,l,theta_2,phi)).real
for m in range(1,l+1):
	fcolors_2 += ((complex(Re_list[m],Im_list[m])*sph_harm(m,l,theta_2,phi)).real) + ((complex(Re_list[m],Im_list[m]).conjugate()*pow(-1,m)*sph_harm(-m,l,theta_2,phi)).real) 	

fmax, fmin = fcolors.max(), fcolors.min()
#fcolors = (fcolors)
#fcolors = fcolors/(fmax-fmin)
fcolors = (fcolors - fmin)/(fmax - fmin)

fmax_2, fmin_2 = fcolors_2.max(), fcolors_2.min()
#fcolors = (fcolors)
#fcolors = fcolors/(fmax-fmin)
fcolors_2 = (fcolors_2 - fmin_2)/(fmax_2 - fmin_2)

# Set the aspect ratio to 1 so our sphere looks spherical
#fig = plt.figure(figsize=plt.figaspect(1.))
fig = plt.figure(figsize=(10,5))
ax = fig.add_subplot(1,2,1, projection='3d')
ax.plot_surface(x, y, z,  rstride=1, cstride=1, facecolors=cm.seismic(fcolors)) 
ax.view_init(45,0)
ax.set_axis_off()
ax = fig.add_subplot(1,2,2, projection='3d')
ax.plot_surface(x, y, z, rstride=1, cstride=1, facecolors=cm.seismic(fcolors)) 
ax.view_init(45,180)
ax.set_axis_off()

plt.title('Parameters are: el_not=%s, tau=%s, lambda3=%s, lambda4=%s\nMinimum Hamiltonian value: %s'%(parameters_list[0], parameters_list[1], parameters_list[2], parameters_list[3], parameters_list[4]), fontsize=8)
#plt.savefig("surface_pattern.png")
plt.show()
