from IPython import embed

import lammpsbuilder as lb
import random
import math


def randomSphere(r):
	x = random.uniform(-1,1)
	y = random.uniform(-1,1)
	z = random.uniform(-1,1)

	scale = r/(math.sqrt(x*x + y*y + z*z))

	x*=scale
	y*=scale
	z*=scale

	return round(x,3),round(y,3),round(z,3)

def voxelise(dim,vox_size):
	for i in range(len(dim)):
		if dim[i] >= 0:
			dim[i] = math.floor(dim[i]/vox_size)*vox_size
		else:
			dim[i] = math.ceil(dim[i]/vox_size)*vox_size
	return dim


def collagenMonomer(length,xlo,xhi,ylo,yhi,zlo,zhi,data,mononomerId,cutoff):

	x_0 = random.uniform(xlo,xhi)
	y_0 = random.uniform(ylo,yhi)
	z_0 = random.uniform(zlo,zhi)

	delta = randomSphere(cutoff)

	charge_dist = [9.994870, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000, 
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000,
				   -9.994750, 9.993550, 0.000000, 0.000000, -0.005370, 0.000000,
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000,
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, -9.999900]

	for i in range(length):
		# occupied_xyz = []
		# for atom in data.atoms:
		# 	occupied_xyz.append([atom.x,atom.y,atom.z])

		# if [x_0+i*delta[0],y_0+i*delta[1],z_0+i*delta[2]] in occupied_xyz:
		#  	print "Duplicate in xyz!", x_0+i*delta[0],  y_0+i*delta[1], z_0+i*delta[2]

		data.addAtom(1,charge_dist[i],x_0+i*delta[0],y_0+i*delta[1],z_0+i*delta[2],mononomerId)
		# embed()
		# data.atoms = data.atoms[:-i+1]

	return 0



def redistribution(data,length):
	ends = []
	for atom in data.atoms:	
		if atom.atomId%length == 1 or atom.atomId%length == 0:
			ends.append(atom)




def main():
	
	xlo = -50
	xhi = 50
	ylo = -50
	yhi = 50
	zlo = -50
	zhi = 50

	excl_zone = 1.12 #should be the lj cut off

	monomer_length = 30

	data = lb.LammpsData(atomTypes=1, xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi, zlo = zlo, zhi = zhi)
	data.addMass(1,1.0)

	#resize box which atoms can be placed in to a multiple of the voxel size (note: does not change actual box size of simulation)
	#dim = voxelise(dim,excl_zone) 

	for i in range(100):
		# collagenMonomer(30,xlo,xhi,ylo,yhi,zlo,zhi,data,i,excl_zone)
		collagenMonomer(monomer_length,xlo,xhi,ylo,yhi,zlo,zhi,data,i,excl_zone)

	redistribution(data, monomer_length)

	f = open('lammps.data', 'w')

	f.write(str(data))

	f.close()

	return 0

if __name__ == "__main__":
    main()