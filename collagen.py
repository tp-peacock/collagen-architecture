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


class Monomer:

	def __init__(self, num_of_atoms, atom_radius, charge_dist):

		self.num_of_atoms = num_of_atoms
		self.atom_radius = atom_radius
		self.length = num_of_atoms*atom_radius*2
		self.charge_dist = charge_dist
		self.end0 = []
		self.end1 = []
		self.delta = None

	def setEndPoints(self,box):

		x0 = random.uniform(box.xlo,box.xhi)
		y0 = random.uniform(box.ylo,box.yhi)
		z0 = random.uniform(box.zlo,box.zhi)
		self.end0 = [x0,y0,z0]

		self.delta = randomSphere(self.atom_radius)
		self.end1 = [x0+self.delta[0]*self.length,y0+self.delta[1]*self.length,z0+self.delta[2]*self.length]
		return 0


class MonomerData(lb.LammpsData):

	def __init__(self,atomTypes=0,bondTypes=0,angleTypes=0,xlo=-200,xhi=200,ylo=-50,yhi=50,zlo=-50,zhi=50):
		lb.LammpsData.__init__(self,atomTypes,bondTypes,angleTypes,xlo,xhi,ylo,yhi,zlo,zhi)
		self.monomers = []


	def addMonomer(self,monomer,monomerId):
		
		for i in range(monomer.num_of_atoms):
			self.addAtom(1,
				         monomer.charge_dist[i],
				         monomer.end0[0]+i*monomer.delta[0],
						 monomer.end0[1]+i*monomer.delta[1],
						 monomer.end0[2]+i*monomer.delta[2],
						 monomerId)
		self.monomers.append(monomer)
		return monomerId


def randomSphere(r):
	x = random.uniform(-1,1)
	y = random.uniform(-1,1)
	z = random.uniform(-1,1)

	scale = r/(math.sqrt(x*x + y*y + z*z))

	x*=scale
	y*=scale
	z*=scale

	return round(x,3),round(y,3),round(z,3)


# def findEnds(data, length):
# 	ends = []
# 	for atom in data.atoms:	
# 		if atom.atomId%length == 1 or atom.atomId%length == 0:
# 			ends.append(atom)	
# 	return ends


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

def closestApproach(monomer1,monomer2):
	print "P0: ", monomer1.end0
	print "P1: ", monomer1.end1
	print "Q0: ", monomer2.end0
	print "Q1: ", monomer2.end1




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


# def redistribute(data, box, atom0,cutoff,monomer_length):
# #	noise = random.randint()*cutoff
# 	idcount = data.atoms[atom0].atomId
# 	print "before: ", data.atoms[atom0]
# 	moleculeId = data.atoms[atom0].moleculeId
# 	for j in range(atom0, atom0+monomer_length):
# 		# data.atoms[i].x -= noise
# 		# data.atoms[i].y -= noise
# 		# data.atoms[i].z -= noise
# 		del data.atoms[atom0]

# 	collagenMonomer(monomer_length,box,data,moleculeId,cutoff)

# 	for k in range(len(data.atoms)-monomer_length,len(data.atoms)):
# 		data.atoms[k].atomId = idcount
# 		idcount+=1

# 	# print "after : ", data.atoms[2970]
# 	print "redistributed"

	#embed()


def main():
	
	box = Box()

	excl_zone = 1.12 #should be the lj cut off

	monomer_length = 30

	data = MonomerData(atomTypes=1, xlo = box.xlo, xhi = box.xhi, ylo = box.ylo, yhi = box.yhi, zlo = box.zlo, zhi = box.zhi)
	data.addMass(1,1.0) 

	# atom1 = ends[0]
	# atom2 = ends[1]
	# close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	# ends = findEnds(data, monomer_length)

	# atom1 = ends[0]
	# atom2 = ends[1]
	# close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	# ends = findEnds(data, monomer_length)
	# atom1 = ends[0]
	# atom2 = ends[1]
	# close_atoms = findCloseAtoms(atom1,atom2,ends,data,excl_zone,monomer_length,box)

	charge_dist = [9.994870, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000, 
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000,
				   -9.994750, 9.993550, 0.000000, 0.000000, -0.005370, 0.000000,
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, 0.000000,
				   0.000000, -0.005370, 0.000000, 0.000000, -0.005370, -9.999900]



	for i in range(20):
		monomer = Monomer(monomer_length, excl_zone,charge_dist)
		monomer.setEndPoints(box)
		data.addMonomer(monomer,i+1)


	closestApproach(data.monomers[0],data.monomers[1])
	

	f = open('lammps.data', 'w')
	f.write(str(data))
	f.close()

	return 0

if __name__ == "__main__":
    main()