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

	return x,y,z



def pythagoras(dim1):
	dim2 = dim1/(math.sqrt(2))
	dim3 = dim1/(math.sqrt(3))
	return dim1, dim2, dim3

def collagenMonomer(length,xlo,xhi,ylo,yhi,zlo,zhi,data,mononomerId):

	x_0 = random.randrange(xlo,xhi,1)
	y_0 = random.randrange(xlo,xhi,1)
	z_0 = random.randrange(xlo,xhi,1)

	# diameters = pythagoras(1)

	# orient = random.randrange(7)

	# for i in range(length):
	# 	delta1d = i*diameters[0]
	# 	delta2d = i*diameters[1]
	# 	delta3d = i*diameters[2]

	# 	if orient == 0:
	# 		data.addAtom(1,0.5,x_0+delta1d,y_0,z_0,mononomerId)
	# 	if orient == 1:
	# 		data.addAtom(1,0.5,x_0,y_0+delta1d,z_0,mononomerId)
	# 	if orient == 2:	
	# 		data.addAtom(1,0.5,x_0,y_0,z_0+delta1d,mononomerId)

	# 	if orient == 3:
	# 		data.addAtom(1,0.5,x_0+delta2d,y_0+delta2d,z_0,mononomerId)
	# 	if orient == 4:
	# 		data.addAtom(1,0.5,x_0,y_0+delta2d,z_0+delta2d,mononomerId)
	# 	if orient == 5:	
	# 		data.addAtom(1,0.5,x_0+delta2d,y_0,z_0+delta2d,mononomerId)

	# 	if orient == 6:	
	# 		data.addAtom(1,0.5,x_0+delta3d,y_0+delta3d,z_0+delta3d,mononomerId)

	delta = randomSphere(1)

	for i in range(length):

		data.addAtom(1,0.5,x_0+i*delta[0],y_0+i*delta[1],z_0+i*delta[2],mononomerId)

	return 0



def main():
	
	xlo = -200
	xhi = 200
	ylo = -50
	yhi = 50
	zlo = -50
	zhi = 50

	data = lb.LammpsData(atomTypes=1, xlo = xlo, xhi = xhi, ylo = ylo, yhi = yhi, zlo = zlo, zhi = zhi)
	data.addMass(1,1.0)

	for i in range(20):	
		collagenMonomer(16,xlo,xhi,ylo,yhi,zlo,zhi,data,i)

	f = open('lammps.data', 'w')

	f.write(str(data))

	f.close()

	for i in range(50):
		print randomSphere(1)

	return 0

if __name__ == "__main__":
    main()