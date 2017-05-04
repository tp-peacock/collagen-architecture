# minimally adapted from lammpsbuilder.py by J C Forster
# https://github.com/takenbymood/membrane-deformation


import sys
import os


class LammpsAtom:

	def __init__(self,atomId,moleculeId,atomType,charge,x=0,y=0,z=0):
		self.atomId = atomId
		self.moleculeId = moleculeId
		self.atomType = atomType
		self.charge = charge
		self.x = x
		self.y = y
		self.z = z

	def __str__(self):
		s=str(self.atomId)+" "+str(self.moleculeId)+" "+str(self.atomType)
		s+=" "+str(self.charge)
		s+=" "+str(self.x)+" "+str(self.y)+" "+str(self.z)
		return s

class LammpsBond:

	def __init__(self,bondId,bondType,atom1,atom2):
		self.bondId = bondId
		self.bondType = bondType
		self.atom1 = atom1
		self.atom2 = atom2

	def __str__(self):
		s=str(self.bondId)+" "+str(self.bondType)+" "+str(self.atom1)+" "+str(self.atom2)
		return(s)

class LammpsAngle:

	def __init__(self,angleId,angleType,atom1,atom2,atom3):
		self.angleId = angleId
		self.angleType = angleType
		self.atom1 = atom1
		self.atom2 = atom2
		self.atom3 = atom3

	def __str__(self):
		s = str(self.angleId)+" "+str(self.angleType)+" "+str(self.atom1)+" "+str(self.atom2)+" "+str(self.atom3)
		return s

class LammpsMass:
	def __init__(self,atomType,mass):
		self.atomType = atomType
		self.mass = mass

	def __str__(self):
		s = str(self.atomType)+" "+str(self.mass)
		return s

class LammpsData:


	def __init__(self,atomTypes=0,bondTypes=0,angleTypes=0,xlo=-200,xhi=200,ylo=-50,yhi=50,zlo=-50,zhi=50):
		self.atomTypes = atomTypes
		self.bondTypes = bondTypes
		self.angleTypes = angleTypes
		self.xlo = xlo
		self.xhi = xhi
		self.ylo = ylo
		self.yhi = yhi
		self.zlo = zlo
		self.zhi = zhi
		self.atoms = []
		self.bonds = []
		self.masses = []
		self.angles = []

	def addAtom(self,atomType,charge,x,y,z=0,moleculeId=-1):
		atomId = len(self.atoms)+1
		if(moleculeId == -1):
			moleculeId = atomId
		a = LammpsAtom(atomId,moleculeId,atomType,charge,x,y,z)
		self.atoms.append(a)
		return atomId

	def addBond(self,bondType,atom1,atom2):
		bondId = len(self.bonds)+1
		bond = LammpsBond(bondId,bondType,atom1,atom2)
		self.bonds.append(bond)
		return bondId

	def addAngle(self,angleType,atom1,atom2,atom3):
		angleId = len(self.angles)+1
		angle = LammpsAngle(angleId,angleType,atom1,atom2,atom3)
		self.angles.append(angle)
		return angleId

	def addMass(self,atomType,mass):
		self.masses.append(LammpsMass(atomType,mass))

	def __str__(self):
		s = "LAMMPS CONFIG FILE\n\n"

		s += "  "+str(len(self.atoms))+" atoms\n"
		s += "  "+str(len(self.bonds))+" bonds\n"
		s += "  "+str(len(self.angles))+" angles\n"
		s += "\n"
		s += "  "+str(self.atomTypes)+" atom types\n"
		s += "  "+str(self.bondTypes)+" bond types\n"
		s += "  "+str(self.angleTypes)+" angle types\n"
		s += "\n"
		s += "  "+str(self.xlo)+"  "+str(self.xhi)+" xlo xhi\n"
		s += "  "+str(self.ylo)+"  "+str(self.yhi)+" ylo yhi\n"
		s += "  "+str(self.zlo)+"  "+str(self.zhi)+" zlo zhi\n"
		s += "\n"
		if(len(self.masses)>0):
			s+="Masses\n\n"
			for m in self.masses:
				s+= "  "+str(m)+"\n"
			s+="\n"
		if(len(self.atoms)>0):
			s+="Atoms\n\n"
			for a in self.atoms:
				s+="  "+str(a)+"\n"
			s+="\n"
		if(len(self.bonds)>0):
			s+="Bonds\n\n"
			for b in self.bonds:
				s+="  "+str(b)+"\n"
			s+="\n"
		if(len(self.angles)>0):
			s+="Angles\n\n"
			for a in self.angles:
				s+="  "+str(a)+"\n"
			s+="\n"
		return s


class LammpsScript:
	sep = "			"

	def addPair(self,atom1,atom2,eps=0,sig=0,cutoff=""):
		s = str(atom1)+" "+str(atom2)+" "+str(eps)+" "+str(sig)+" "+str(cutoff)
		self.pair_coeffs.append(s)

	def addPairModify(self,line):
		self.pairmod.append(line)

	def addBond(self,bond,K,x0):
		s = str(bond)+" "+str(K)+" "+str(x0)
		self.bond_coeffs.append(s)

	def addAngle(self,angle,K,theta0):
		s = str(angle)+" "+str(K)+" "+str(theta0)
		self.angle_coeffs.append(s)

	def addGroup(self,name,members,order="type"):
		s = str(name)+" "+str(order)+" "
		for m in members:
			s+=str(m)+" "
		self.groups.append(s)

	def addFix(self,group,action):
		s = str(group)+ " " + str(action)
		self.fixes.append(s)

	def addLine(self,line):
		self.lines.append(str(line))

	def __init__(self,read_data="",dump="id all xyz 100 out.xyz",thermo="300",timestep="0.001",run="100000",dimension="2",units="lj",velocity="all create 1.0 1000",atom_style="molecular",atom_modify="sort 0 1",neighbour="0.3 bin",neigh_modify="every 1 delay 1",angle_style="harmonic",bond_style="harmonic",pair_style="lj/cut 2.5"):
		self.read_data = read_data
		self.dump = dump
		self.dimension = dimension
		self.units = units
		self.atom_style = atom_style
		self.atom_modify = atom_modify
		self.neighbour = neighbour
		self.neigh_modify = neigh_modify
		self.angle_style = angle_style
		self.bond_style = bond_style
		self.pair_style = pair_style
		self.timestep = timestep
		self.thermo = thermo
		self.run = run
		self.velocity = velocity
		self.bond_coeffs=[]
		self.pair_coeffs=[]
		self.angle_coeffs=[]
		self.groups = []
		self.fixes = []
		self.lines = []
		self.pairmod=[]

	def __str__(self):
		d = self.sep
		s="dimension"+d+self.dimension+"\n"
		s+="units"+d+self.units+"\n"
		s+="atom_style"+d+self.atom_style+"\n"
		s+="boundary f f p\n"
		s+="atom_modify"+d+self.atom_modify+"\n"
		s+="\n"
		s+="read_data"+d+self.read_data+"\n"
		#s+="neighbour"+d+self.neighbour+"\n"
		s+="neigh_modify"+d+self.neigh_modify+"\n"
		if len(self.bond_coeffs)>0:
			s+="bond_style"+d+self.bond_style+"\n"
		for b in self.bond_coeffs:
			s+="bond_coeff"+d+b+"\n"
		if len(self.angle_coeffs)>0:
			s+="angle_style"+d+self.angle_style+"\n"
		for a in self.angle_coeffs:
			s+="angle_coeff"+d+a+"\n"
		if len(self.pair_coeffs)>0:
			s+="pair_style"+d+self.pair_style+"\n"
		for p in self.pair_coeffs:
			s+="pair_coeff"+d+p+"\n"
		for p in self.pairmod:
			s+="pair_modify"+d+p+"\n"
		
		s+="dump"+d+self.dump+"\n"

		for g in self.groups:
			s+="group"+d+g+"\n"

		s+="velocity move create 1.0 1\n"
		s+="velocity protein set 0 -4 0\n"
		i=0
		for f in self.fixes:
			s+="fix"+d+str(i)+" "+f+"\n"
			i+=1
		s+="\n"
		for l in self.lines:
			s+=l
		s+="\n"
		s+="thermo"+d+self.thermo+"\n"
		s+="timestep"+d+self.timestep+"\n"
		s+="run"+d+self.run+"\n"

		return s
		

class LammpsSimulation:

	def __init__(self,name,filedir="",run="100000"):

		self.name = name
		self.scriptName = name+"_script.in"
		self.dataName = name+"_data.data"
		self.filedir = filedir
		self.script = LammpsScript(read_data=self.dataName,run=run)
		self.data = LammpsData()

	def __del__(self):
		name = ""
		scriptName = ""
		dataName = ""
		filedir = ""
		script = None
		data = None

	def saveFiles(self):
		with open(self.filedir+self.scriptName, 'w') as file_:
			file_.write(str(self.script))

		with open(self.filedir+self.dataName, 'w') as file_:
			file_.write(str(self.data))

	def deleteFiles(self):
		os.remove(self.filedir+self.scriptName)
		os.remove(self.filedir+self.dataName)

	def __str__(self):
		return self.name + "\n" + str(self.script) + "\n" + str(self.data)