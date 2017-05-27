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

	def copy(self, monomer):
		self.__init__(self,monomer.num_of_atoms,monomer.atom_radius,monomer.charge_dist)

	def setEndPoints(self,box):

		x0 = random.uniform(box.xlo, box.xhi)
		y0 = random.uniform(box.ylo, box.yhi)
		z0 = random.uniform(box.zlo, box.zhi)
		self.end0 = [x0,y0,z0]

		self.delta = randomSphere(self.atom_radius)
		#self.end1 = [None,None,None]

		# if x0 < (box.xhi + box.xlo)/2:
		# 	self.end1[0] = x0+abs(self.delta[0]*self.length)
		# else:
		# 	self.end1[0] = x0-abs(self.delta[0]*self.length)

		# if y0 < (box.yhi + box.ylo)/2:
		# 	self.end1[1] = y0+abs(self.delta[1]*self.length)
		# else:
		# 	self.end1[1] = y0-abs(self.delta[1]*self.length)

		# if z0 < (box.zhi + box.zlo)/2:
		# 	self.end1[2] = z0+abs(self.delta[2]*self.length)
		# else:
		# 	self.end1[2] = z0-abs(self.delta[2]*self.length)

		
		self.end1 = [x0+self.delta[0]*self.length,y0+self.delta[1]*self.length,z0+self.delta[2]*self.length]

		if self.end1[0] > box.xhi - 2.0:
		#	embed()
			self.end0[0] -= self.end1[0] - (box.xhi - 2.0)
			self.end1[0] = box.xhi-2.0

		if self.end1[0] < box.xlo + 2.0:
			self.end0[0] -= self.end1[0] - (box.xlo + 2.0)
			self.end1[0] = box.xlo +2.0

		if self.end1[1] > box.yhi - 2.0:
			self.end0[1] -= self.end1[1] - (box.yhi - 2.0)
			self.end1[1] = box.yhi - 2.0

		if self.end1[1] < box.ylo + 2.0:
			self.end0[1] -= self.end1[1] - (box.ylo + 2.0)
			self.end1[1] = box.ylo + 2.0

		if self.end1[2] > box.zhi - 2.0:
			self.end0[2] -= self.end1[2] - (box.zhi - 2.0)
			self.end1[2] = box.zhi - 2.0

		if self.end1[2] < box.zlo + 2.0:
			self.end0[2] -= self.end1[2] - (box.zlo + 2.0)
			self.end1[2] = box.zlo + 2.0


		return 0
		# print "end0: ", self.end0
		# print "end1: ", self.end1
		# print "x_diff", self.end1[0]  -self.end0[0]
		# print "y_diff", self.end1[1]  -self.end0[1]
		# print "z_diff", self.end1[2]  -self.end0[2]
		# print "length", np.linalg.norm(np.array(self.end1)-np.array(self.end0))
		# print self.length

		
		# for i in range(len(self.end1)):
		# 	if self.end1[i] > 50.0 or self.end1[i] < -50.0:
		# 		embed()

	#	return 0


class MonomerData(lb.LammpsData):

	def __init__(self,atomTypes=0,bondTypes=0,angleTypes=0,xlo=-200,xhi=200,ylo=-50,yhi=50,zlo=-50,zhi=50):
		lb.LammpsData.__init__(self,atomTypes,bondTypes,angleTypes,xlo,xhi,ylo,yhi,zlo,zhi)
		self.monomers = []


	def addMonomer(self,monomer,monomerId):
		
		for i in range(monomer.num_of_atoms):

			# if monomer.end1[0] >= monomer.end0[0]:	
			# 	proj_x = monomer.end0[0]+abs(i*monomer.delta[0])
			# else:
			# 	proj_x = monomer.end0[0]-abs(i*monomer.delta[0])

			# if monomer.end1[1] >= monomer.end0[1]:
			# 	proj_y = monomer.end0[1]+abs(i*monomer.delta[1])
			# else:
			# 	proj_y = monomer.end0[1]-abs(i*monomer.delta[1])

			# if monomer.end1[2] >= monomer.end0[2]:
			# 	proj_z = monomer.end0[2]+abs(i*monomer.delta[2])
			# else:
			# 	proj_z = monomer.end0[2]-abs(i*monomer.delta[2])

			# x = proj_x
			# y = proj_y
			# z = proj_z

			x = monomer.end0[0]+i*monomer.delta[0]	
			y = monomer.end0[1]+i*monomer.delta[1]
			z = monomer.end0[2]+i*monomer.delta[2]		

			self.addAtom(1,
				         monomer.charge_dist[i],
				         x,
						 y,
						 z,
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

def clamp(A0,A1,Ac):
	for i in range(3):
		if (Ac[i] < A0[i] and Ac[i] < A1[i]) or (Ac[i] > A0[i] and Ac[i] > A1[i]):
			if min(abs(A1[i]-Ac[i]),abs(A0[i]-Ac[i])) == abs(A1[i]-Ac[i]):
				Ac[i] = A1[i]
			else:
				Ac[i] = A0[i]
	return Ac

def closestApproach(monomer1,monomer2):
	P0 = np.array(monomer1.end0)
	P1 = np.array(monomer1.end1)
	Q0 = np.array(monomer2.end0)
	Q1 = np.array(monomer2.end1)

	u = (P1 - P0)/np.linalg.norm(P1-P0)
	v = (Q1 - Q0)/np.linalg.norm(Q1-Q0)

	a = u.dot(u)
	b = u.dot(v)
	c = v.dot(v)
	d = u.dot(P0-Q0)
	e = v.dot(P0-Q0)

	Sc = (b*e - c*d)/(a*c - b*b)
	Tc = (a*e - b*d)/(a*c - b*b)

	Pc = P0 + Sc*u
	Qc = Q0 + Tc*v

	Pc = clamp(P0,P1,Pc)
	Qc = clamp(Q0,Q1,Qc)

	CA = np.linalg.norm(Pc - Qc)

	return CA

def isWithinProximity(monomer1,monomer2,proximity):
	return True if closestApproach(monomer1,monomer2) <= proximity else False

def redistribute(monomer, noiselimit):
	noise = random.uniform(0.5*noiselimit,1.5*noiselimit)
	sign = -1 if random.randint(0,1) == 0 else 1
	for i in range(3):
		monomer.end0[i]+=sign*noise
		monomer.end1[i]+=sign*noise

def buildNewMonomer(monomer,box):
	new_mon = Monomer(monomer.num_of_atoms, monomer.atom_radius, monomer.charge_dist)
	new_mon.setEndPoints(box)
	return new_mon

def replace_monomer(empty_monomers, monomer_length, atom_radius, charge_dist, box, offender, excl_zone):
	new_mon = Monomer(monomer_length, atom_radius, charge_dist)
	new_mon.setEndPoints(box)
	empty_monomers[offender] = new_mon
	for i in range(offender):
		CA = closestApproach(empty_monomers[i],empty_monomers[offender])
		if CA <= 3*excl_zone:
			print "Overlap!"
			replace_monomer(empty_monomers, monomer_length, atom_radius, charge_dist, box, offender, excl_zone)


def overlapsWithAny(monomer1,monomerIndex,monomers,proximity):
	for i in range(len(monomers)):
		if i!=monomerIndex:
			monomer2 = monomers[i]
			if isWithinProximity(monomer1,monomer2,proximity):
				return True
	return False

def fixOverlaps(monomer1,monomerIndex,monomers,proximity,box):
	if not overlapsWithAny(monomer1,monomerIndex,monomers,proximity) :
		return monomers
	else:
		monomers[monomerIndex] = buildNewMonomer(monomer1,box)
		return fixOverlaps(monomers[monomerIndex],monomerIndex,monomers,proximity,box)


def main():
	
	box = Box()

	num_of_monomers = 100
	excl_zone = 1.3 #should be the lj cut off
	atom_radius = 0.5

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

	redistribution_noise = 1
	redist_counter = 0
	empty_monomers = []

	for i in range(num_of_monomers):
		monomer = Monomer(monomer_length, atom_radius,charge_dist)
		monomer.setEndPoints(box)
	#	data.addMonomer(monomer,i+1)
		empty_monomers.append(monomer)

	overlap = None
	offender = None

	for i in range(num_of_monomers):
		print "building: ", i
		monomer1 = empty_monomers[i]
		empty_monomers = fixOverlaps(monomer1,i,empty_monomers,2*excl_zone,box)


	# while overlap != 0:
	# 	overlap = 0
	# 	for i in range(num_of_monomers):
	# 		for j in range(num_of_monomers):
	# 			if i < j
	# 			#	CA = closestApproach(data.monomers[i],data.monomers[j])
	# 				CA = closestApproach(empty_monomers[i],empty_monomers[j])
	# 				if CA <= 3*excl_zone:
	# 					overlap += 1
	# 					offender = j
	# 					print "Overlap!"
	# 					break
	# 		if overlap > 0:
	# 			break

	# 	if overlap > 0:
	# 		redist_counter+=1
	# 		print "comparing i: ", i, " offender: ", offender
	# 		# empty_monomers.remove(empty_monomers[j])
	# 		# new_mon = Monomer(monomer_length, atom_radius,charge_dist)
	# 		# new_mon.setEndPoints(box)
	# 		# empty_monomers.append(new_mon)

	# 		replace_monomer(empty_monomers, monomer_length, atom_radius, charge_dist, box, offender, excl_zone)


	# 	#	redistribute(data.monomers[offender],redistribution_noise)

	# 		overlap = None



	for i in range(len(empty_monomers)):
		data.addMonomer(empty_monomers[i],i+1)
	

#	f = open('lammps.data', 'w')
	f = open('small_sim.data', 'w')
	f.write(str(data))
	f.close()

	return 0

if __name__ == "__main__":
    main()