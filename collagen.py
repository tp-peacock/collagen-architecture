from IPython import embed

import lammpsbuilder as lb
import random
import math
import numpy as np

class Box:

	def __init__(self,xlo=-50,xhi=50,ylo=-50,yhi=50,zlo=-50,zhi=50):

		self.xlo = xlo
		self.xhi = xhi
		self.ylo = ylo
		self.yhi = yhi
		self.zlo = zlo
		self.zhi = zhi

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


def collagenMonomer(length,box,data,mononomerId,cutoff):

	x_0 = random.uniform(box.xlo,box.xhi)
	y_0 = random.uniform(box.ylo,box.yhi)
	z_0 = random.uniform(box.zlo,box.zhi)

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


def findEnds(data, length):
	ends = []
	for atom in data.atoms:	
		if atom.atomId%length == 1 or atom.atomId%length == 0:
			ends.append(atom)	
	return ends


def closestApproach(ends):

	P0 = [ends[0].x,ends[0].y,ends[0].z]
	Pe = [ends[1].x,ends[1].y,ends[1].z]
	Q0 = [ends[2].x,ends[2].y,ends[2].z]
	Qe = [ends[3].x,ends[3].y,ends[3].z]

	u = []
	v = []
	xp = Pe[0]-P0[0]
	yp = Pe[1]-P0[1]
	zp = Pe[2]-P0[2]
	xq = Qe[0]-Q0[0]
	yq = Qe[1]-Q0[1]
	zq = Qe[2]-Q0[2]

	for i in range(3):
		u.append((Pe[i]-P0[i])/math.sqrt( xp*xp + yp*yp + zp*zp ))

	for i in range(3):
		v.append((Qe[i]-Q0[i])/math.sqrt( xq*xq + yq*yq + zq*zq ))

	a = u[0]*u[0] + u[1]*u[1] + u[2]*u[2]
	b = u[0]*v[0] + u[1]*v[1] + u[2]*v[2]
	c = v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
	d = u[0]*(P0[0]-Q0[0]) + u[1]*(P0[1]-Q0[1]) + u[2]*(P0[2]-Q0[2])
	e = v[0]*(P0[0]-Q0[0]) + v[1]*(P0[1]-Q0[1]) + v[2]*(P0[2]-Q0[2])

	u_np = np.array(u)
	v_np = np.array(v)
	P0_np = np.array(P0)
	Q0_np = np.array(Q0)

	closest = (P0_np - Q0_np) + ((b*e - c*d)*u_np - (a*e -b*d)*v_np)/(a*c - b*b)

	abs_dist = np.linalg.norm(closest)

def cApproach(atom1, atom2, atom3, atom4):

	P0 = np.array([atom1.x, atom1.y, atom1.z])
	Pe = np.array([atom2.x, atom2.y, atom2.z])
	Q0 = np.array([atom3.x, atom3.y, atom3.z])
	Qe = np.array([atom4.x, atom4.y, atom4.z])

	u = (Pe - P0)/np.linalg.norm(Pe-P0)
	v = (Qe - Q0)/np.linalg.norm(Qe-Q0)

	a = u.dot(u)
	b = u.dot(v)
	c = v.dot(v)
	d = u.dot(P0-Q0)
	e = v.dot(P0-Q0)

	closest = np.linalg.norm( (P0 - Q0) + ((b*e - c*d)*u - (a*e -b*d)*v)/(a*c - b*b) )

	return closest


def findCloseAtoms(atom1, atom2,ends,data,cutoff,monomer_length,box):

	close_atoms = []

	flag = False
	while flag == False:
		for i in range(2, len(ends),2):
			
			flag = True
			closest = cApproach(atom1, atom2, ends[i], ends[i+1])
		
			if closest <= cutoff*2:
				close_atoms.append(i*monomer_length/2)
				redistribute(data, box, i*monomer_length/2,cutoff,monomer_length)
				ends = findEnds(data,monomer_length)
				flag = False
				break

	print "close atoms: ", close_atoms
	return close_atoms


def redistribute(data, box, atom0,cutoff,monomer_length):
#	noise = random.randint()*cutoff
	idcount = data.atoms[atom0].atomId
	print "before: ", data.atoms[atom0]
	moleculeId = data.atoms[atom0].moleculeId
	for j in range(atom0, atom0+monomer_length):
		# data.atoms[i].x -= noise
		# data.atoms[i].y -= noise
		# data.atoms[i].z -= noise
		del data.atoms[atom0]

	collagenMonomer(monomer_length,box,data,moleculeId,cutoff)

	for k in range(len(data.atoms)-monomer_length,len(data.atoms)):
		data.atoms[k].atomId = idcount
		idcount+=1

	print "after : ", data.atoms[2970]
	print "redistributed"

	#embed()


def main():
	
	box = Box()

	excl_zone = 1.12 #should be the lj cut off

	monomer_length = 30

	data = lb.LammpsData(atomTypes=1, xlo = box.xlo, xhi = box.xhi, ylo = box.ylo, yhi = box.yhi, zlo = box.zlo, zhi = box.zhi)
	data.addMass(1,1.0)

	#resize box which atoms can be placed in to a multiple of the voxel size (note: does not change actual box size of simulation)
	#dim = voxelise(dim,excl_zone) 

	for i in range(100):
		# collagenMonomer(30,xlo,xhi,ylo,yhi,zlo,zhi,data,i,excl_zone)
		collagenMonomer(monomer_length,box,data,i,excl_zone)

	ends = findEnds(data, monomer_length)

	atom1 = ends[0]
	atom2 = ends[1]
	close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	ends = findEnds(data, monomer_length)

	atom1 = ends[0]
	atom2 = ends[1]
	close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	ends = findEnds(data, monomer_length)
	atom1 = ends[0]
	atom2 = ends[1]
	close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	#redistribute(close_atoms)




	f = open('lammps.data', 'w')

	f.write(str(data))

	f.close()

	return 0

if __name__ == "__main__":
    main()