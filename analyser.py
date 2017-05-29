import sys
from IPython import embed
import numpy as np

class Monomer:

	def __init__(self, index, num_of_atoms, atom_radius):
		self.num_of_atoms = num_of_atoms
		self.atom_radius = atom_radius
		self.length = num_of_atoms*atom_radius*2
		self.end0 = []
		self.end1 = []
		self.index = index


def getFinalState(filename):

	f = open(filename, 'r')
	
	lines = []
	lcount = 0
	for line in f:
		lines.append(line)
		if "Timestep" in line:
			finalstep = lcount
		lcount += 1
	f.close()
	lines = lines[finalstep+1:]

	return lines

def buildFinalStateFile(filename, lines):

	f = open(filename, 'w')
	f.write(str(len(lines))+"\n")
	f.write("Atoms. Timestep: 0\n")
	print "Writing final state file..."
	for i in range(len(lines)):
		f.write(lines[i])
	f.close()
	print "Final state file written to", filename
	return 0




def buildMonomers(lines,monomer_num_of_atoms, atom_radius):
	data = []
	index = 0
	for i in range(len(lines)):
		if i % 30 == 0:
			new_mon = Monomer(index,monomer_num_of_atoms,atom_radius)

			new_mon.end0 = [float(lines[i].split(" ")[1]),float(lines[i].split(" ")[2]),float(lines[i].split(" ")[3][:-1])]
			new_mon.end1 = [float(lines[i+monomer_num_of_atoms-1].split(" ")[1]),float(lines[i+monomer_num_of_atoms-1].split(" ")[2]),float(lines[i+monomer_num_of_atoms-1].split(" ")[3][:-1])]
			data.append(new_mon)
			index+=1
	return(data)

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

def findNeighbours(data,i):

	neighbours = [data[i]]
	for j in range(len(data)):
		if j == i:
			continue
		else:
			CA = closestApproach(data[i], data[j])
			if CA <= 1:
				neighbours.append(data[j])

	return neighbours

def buildNeighbourFile(lines,data,filename,n, monomer_num_of_atoms):
	neighbours = findNeighbours(data,n)
	f = open(filename, 'w')
	f.write(str(len(neighbours)*monomer_num_of_atoms)+'\n')
	f.write('Atoms. Timestep: 0\n')
	for i in range(len(neighbours)):
		start_atom = (neighbours[i].index)*monomer_num_of_atoms
		for j in range(monomer_num_of_atoms):
			f.write(lines[start_atom+j])
	f.close()


def main():

	monomer_num_of_atoms = 30
	atom_radius = 0.5

	filename = sys.argv[1]

	lines = getFinalState(filename)

	if '-fs' in sys.argv:
		buildFinalStateFile("final_state.xyz", lines)

	data = buildMonomers(lines, monomer_num_of_atoms, atom_radius)

	if '-n' in sys.argv:
		n = sys.argv.index('-n')
		buildNeighbourFile(lines,data,"neighbours_of_"+str(sys.argv[n-1])+".xyz", int(sys.argv[n-1]), monomer_num_of_atoms)

	findNeighbours(data,0)



	return 0

if __name__ == "__main__":
    main()