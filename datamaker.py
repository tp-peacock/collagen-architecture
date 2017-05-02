"""
-------------------------------------------------------------------------------------------------------
-------------------------------------------------------------------------------------------------------
Changeable Initial Conditions
-------------------------------------------------------------------------------------------------------
"""

# Options - specify name, number of component, and number of types of component
Atoms = ["atoms", 20, 1]
Bonds = ["bonds", 0, 0]
Angles = ["angles", 0, 0]
Dihedrals = ["dihedrals", 0, 0]
Impropers = ["impropers", 0, 0]

# Provide the masses for each type of atom. Masses array should have the same no. of elements as
# there a
Masses = [1]

#Change the box boundaries
X_bounds = [-10.0,11.0]
Y_bounds = [-12.0,13.0]
Z_bounds = [-14.0,15.0]

"""
-------------------------------------------------------------------------------------------------------
"""

Components = [Atoms, Bonds, Angles, Dihedrals, Impropers]
Box = [X_bounds,Y_bounds,Z_bounds]


from IPython import embed
#embed()

#Initialising options
def initialiser(components):
	num = []
	types = []
	for comp in components:
		if comp[1] != 0:
			num.append(str(comp[1])+" "+comp[0]+"\n")
			types.append(str(comp[2])+" "+comp[0][:-1] +" types\n")
	return [num, types]

def mass_config(masses):
	config = ["Masses\n\n"]
	m = 0
	for mass in masses:
		m += 1
		config.append(str(m) +" "+str(mass)+"\n")
	return config

def box_config(box):
	config = ["xlo xhi", "ylo yhi", "zlo zhi"]
	counter = 0
	for coords in box:
		config[counter] = str(coords[0])+" "+str(coords[1])+" "+config[counter]+"\n"
		counter+=1
	return config

def main():
	#Global Checks
	if len(Masses) != Atoms[2]:
		print "Error! Number of Masses and Atom Types do not match!\n"
		return 0

	f = open('config.data', 'w')

	f.write('LAMMPS Config File\n')
	f.write("\n")
	
	#number of atoms, bonds, angles, dihedrals, impropers and the number of types for each
	comp_options = initialiser(Components)
	for line in comp_options[0]:
		f.write(line)
	f.write ("\n")
	for line in comp_options[1]:
		f.write(line)
	f.write("\n")

	#box configuration
	box_options = box_config(Box)
	for line in box_options:
		f.write(line)
	f.write("\n")

	#mass configuration
	mass_options = mass_config(Masses)
	for line in mass_options:
		f.write(line)
	f.write("\n")

	f.close()

	return 0

if __name__ == "__main__":
    main()