import sys
from IPython import embed
import numpy as np
import matplotlib.pyplot as plt


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
			if CA <= 1.3:
				neighbours.append(data[j])

	return neighbours


def findAlignment(monomer1, monomer2):
	P0 = np.array(monomer1.end0)
	P1 = np.array(monomer1.end1)
	Q0 = np.array(monomer2.end0)
	Q1 = np.array(monomer2.end1)

	# u = (P1 - P0)/np.linalg.norm(P1-P0)
	# v = (Q1 - Q0)/np.linalg.norm(Q1-Q0)

	u = (P1 - P0)
	v = (Q1 - Q0)

	dproduct = u.dot(v)
	mproduct = np.linalg.norm(u)*np.linalg.norm(v)

#	print dproduct, mproduct

	error = 1.5

	if dproduct >= mproduct - error and dproduct <= mproduct + error:
		return "aligned"
	elif dproduct >= -mproduct - error and dproduct <= -mproduct + error:
		return "anti-aligned"
	else:
		return "not aligned"


def findFibres(neighbours):
	aligned = [neighbours[0]]
	align_count = 0
	anti_align_count = 0

	align_d_spacing = []
	anti_align_d_spacing = []

	for i in range(1,len(neighbours)):

		neigh = findAlignment(neighbours[0],neighbours[i]) 
		if neigh == "aligned" or neigh == "anti-aligned":
			aligned.append(neighbours[i])
		if neigh == "aligned":
			align_count += 1
			dspacing = findDSpacing(neighbours[0],neighbours[i])
			align_d_spacing.append(dspacing)
		if neigh == "anti-aligned":
			anti_align_count += 1
			dspacing = findDSpacing(neighbours[0],neighbours[i])
			anti_align_d_spacing.append(dspacing)

	return aligned, align_count, anti_align_count, align_d_spacing, anti_align_d_spacing


def findDSpacing(monomer1, monomer2):
	P0 = np.array(monomer1.end0)
	P1 = np.array(monomer1.end1)
	Q0 = np.array(monomer2.end0)
	Q1 = np.array(monomer2.end1)

	d_vector = (P1 + P0)/2 - (Q1 + Q0)/2
	d = np.linalg.norm(d_vector)

	length = np.linalg.norm(P1 - P0)

	d_spacing = d*30/length

	return d_spacing

def buildNeighbourFile(lines,data,filename,neighbours, monomer_num_of_atoms):
	f = open(filename, 'w')
	f.write(str(len(neighbours)*monomer_num_of_atoms)+'\n')
	f.write('Atoms. Timestep: 0\n')
	for i in range(len(neighbours)):
		start_atom = (neighbours[i].index)*monomer_num_of_atoms
		for j in range(monomer_num_of_atoms):
			f.write(lines[start_atom+j])
	f.close()


def buildPlots(raw_spacings):
	int_spacings=[]
	spacings = [item for sublist in raw_spacings for item in sublist] 
	for i in range(len(spacings)):
		int_spacings.append(int(round(spacings[i])))

	x = []
	y = []
	
	for i in range(1,max(int_spacings)):
		x.append(i)
		y.append(int_spacings.count(i))

	plot = 0



	return plot, int_spacings, x, y



def main():

	monomer_num_of_atoms = 30
	atom_radius = 0.5

	filename = sys.argv[1]

	lines = getFinalState(filename)

	if '-fs' in sys.argv:
		buildFinalStateFile("final_state.xyz", lines)

	data = buildMonomers(lines, monomer_num_of_atoms, atom_radius)




	if '-n' in sys.argv and '-a' in sys.argv:
		n = int(sys.argv[sys.argv.index('-n')-1])
		print "n: ", n
		neighbours = findNeighbours(data,n)

		aligned = findFibres(neighbours)[0]
		buildNeighbourFile(lines,data,"neighbours_of_"+str(n)+".xyz", aligned, monomer_num_of_atoms)


	if '-n' not in sys.argv and '-a' in sys.argv:
		n = int(sys.argv[sys.argv.index('-a')-1])
		print "n: ", n
		neighbours = findNeighbours(data,n)

		aligned = findFibres(neighbours)[0]


	if '-n' in sys.argv and '-a' not in sys.argv:
		n = int(sys.argv[sys.argv.index('-n')-1])
		neighbours = findNeighbours(data,n)
		buildNeighbourFile(lines,data,"neighbours_of_"+str(n)+".xyz", neighbours, monomer_num_of_atoms)


	total_align_count = []
	total_anti_align_count = []
	aligned_d_spacings = []
	anti_aligned_d_spacings = []

	for z in range(0,2000):
	 	neighbours = findNeighbours(data,z)
	 	aligned = findFibres(neighbours)
	 	total_align_count.append(aligned[1])
	 	total_anti_align_count.append(aligned[2])
	 	aligned_d_spacings.append(aligned[3])
	 	anti_aligned_d_spacings.append(aligned[4])
#	 	print aligned[0], aligned[1], aligned[2]


	running_align_total = 0
	running_anti_align_total = 0

#	print max(total_align_count)
#	print max(total_anti_align_count)

	for i in range(max(total_align_count),0,-1):
		running_align_total += total_align_count.count(i)
		print "No. that have ", i, "aligned neighbours", running_align_total

	for i in range(max(total_anti_align_count),0,-1):
		running_anti_align_total += total_anti_align_count.count(i)
		print "No. that have", i, "anti-aligned neighbours", running_anti_align_total



	plot_data1 = buildPlots(aligned_d_spacings)
	x1 = plot_data1[2]
	y1 = plot_data1[3]

	plot_data2 = buildPlots(anti_aligned_d_spacings)
	x2 = plot_data2[2]
	y2 = plot_data2[3]


	fig = plt.figure()


	ax = fig.add_subplot(111)
	ax.plot(x1,y1)
	ax.plot(x2,y2)
	ax.spines['left'].set_position('zero')
	ax.spines['right'].set_color('none')
	ax.spines['bottom'].set_position('zero')
	ax.spines['top'].set_color('none')

# remove the ticks from the top and right edges
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	plt.ylabel('Frequency')
	plt.xlabel('D Spacing (number of atoms)')
	

	ax.set_xlim([0,30])

	plt.savefig('test.png')








#	plot1 = buildPlots(aligned_d_spacings)	
#	plot2 = buildPlots(anti_aligned_d_spacings)

	embed()
	# 	if len(aligned) > 2:
	# 		print "z: ", z, "aligned: ", len(aligned)



	print "number with aligned neighbours:      ", total_align_count
	print "number with anti-aligned neighbours: ", total_anti_align_count


	return 0

if __name__ == "__main__":
    main()