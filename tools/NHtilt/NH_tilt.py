#!/usr/local/bin/python
import sys
from math import *
import numpy
######################################
# ROUTINES TO CALCULATE EULER ANGLES #
######################################
# define the two versors parallel to N-C1 and N-C2 bonds
def peptidePlaneVersors(ncoord,c1coord,c2coord):
	v1x = c1coord[0] - ncoord[0]
	v1y = c1coord[1] - ncoord[1]
	v1z = c1coord[2] - ncoord[2]
	norm = sqrt(v1x*v1x + v1y*v1y + v1z*v1z)
	v1x, v1y, v1z = v1x/norm, v1y/norm, v1z/norm
	v2x = c2coord[0] - ncoord[0]
	v2y = c2coord[1] - ncoord[1]
	v2z = c2coord[2] - ncoord[2]
	norm = sqrt(v2x*v2x + v2y*v2y + v2z*v2z)
	v2x, v2y, v2z = v2x/norm, v2y/norm, v2z/norm
	return (v1x,v1y,v1z,v2x,v2y,v2z)
# calculate e1 and e2 versors
def e1e2 (x3,y3,z3,v1x,v1y,v1z):
	# find e1 = e3 x v1
	x1 = y3*v1z - z3*v1y
	y1 = z3*v1x - x3*v1z
	z1 = x3*v1y - y3*v1x
	norm = sqrt(x1*x1 + y1*y1 + z1*z1)
	x1, y1, z1 = x1/norm, y1/norm, z1/norm
	# find e2 = e3 x e1
	x2 = y3*z1 - z3*y1
	y2 = z3*x1 - x3*z1
	z2 = x3*y1 - y3*x1
	return (x1,y1,z1,x2,y2,z2)
# method 1
def Z_NH(ncoord,c0coord,c1coord,c2coord):
	v1x,v1y,v1z,v2x,v2y,v2z = peptidePlaneVersors(ncoord,c1coord,c2coord)
	# find e3 parallel to N-H bond (which is considered as the bisector of the C-N-Ca angle
	x3 = -(v1x + v2x)
	y3 = -(v1y + v2y)
	z3 = -(v1z + v2z)
	norm = sqrt(x3*x3 + y3*y3 + z3*z3)
	x3, y3, z3 = x3/norm, y3/norm, z3/norm
	# find e1 and e2
	x1,y1,z1,x2,y2,z2 = e1e2(x3,y3,z3,-v1x,-v1y,-v1z)
	#**********************************************
	# DEBUGGING INFO : director cosines per residue
	#**********************************************
	##print nres
	##print x1, y1, z1
	##print x2, y2, z2
	##print x3, y3, z3
	#**********************************************
	R = numpy.array([[-x2,-y2,-z2],[x1,y1,z1],[x3,y3,z3]]) # X and Y are swapped
	return (R,ncoord[0],ncoord[1],ncoord[2])
# method 2
def Z_CC(ncoord,c0coord,c1coord,c2coord):
	v1x,v1y,v1z,v2x,v2y,v2z = peptidePlaneVersors(ncoord,c1coord,c2coord)
	# find e3 parallel to C(a)'-C(a) direction
	x3 = c2coord[0] - c0coord[0]
	y3 = c2coord[1] - c0coord[1]
	z3 = c2coord[2] - c0coord[2]
	norm = sqrt(x3*x3 + y3*y3 + z3*z3)
	x3, y3, z3 = x3/norm, y3/norm, z3/norm
	# find e1 and e2
	x1,y1,z1,x2,y2,z2 = e1e2(x3,y3,z3,v1x,v1y,v1z)
	#**********************************************
	# DEBUGGING INFO : director cosines per residue
	#**********************************************
	##print nres
	##print x1, y1, z1
	##print x2, y2, z2
	##print x3, y3, z3
	#**********************************************
	R = numpy.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
	return (R,c0coord[0],c0coord[1],c0coord[2])
# method 3
def Z_CN(ncoord,c0coord,c1coord,c2coord):
	v1x,v1y,v1z,v2x,v2y,v2z = peptidePlaneVersors(ncoord,c1coord,c2coord)
	# find e3 parallel to C-N bond
	x3, y3, z3 = -v1x, -v1y, -v1z
	# find e1 and e2
	x1,y1,z1,x2,y2,z2 = e1e2(x3,y3,z3,v2x,v2y,v2z)
	#**********************************************
	# DEBUGGING INFO : director cosines per residue
	#**********************************************
	##print nres
	##print x1, y1, z1
	##print x2, y2, z2
	##print x3, y3, z3
	#**********************************************
	R = numpy.array([[x1,y1,z1],[x2,y2,z2],[x3,y3,z3]])
	return (R,c1coord[0],c1coord[1],c1coord[2])
# unknown method error message
def UNKNOWN_METHOD(ncoord,c0coord,c1coord,c2coord):
	print ">>> ERROR <<<"
	print "Method not recognized. Choose among 1, 2 or 3"
	print ""
	exit(1)
#################
# HELLO MESSAGE #
#################
print ""
print "============================================================"
print "=                           N-H TILT                       ="
print "=                                                          ="
print "= Calculates Euler angles describing the orientation of    ="
print "= the potential director frame with respect to the protein ="
print "= diffusion frame. The following choices are available for ="
print "= the definition of the director frame:                    ="
print "=                                                          ="
print "= 1) Z || to N-H bond, Y _|_ peptide plane                 ="
print "=                                                          ="
print "= 2) Z || to C(a)-C(a)' direction, Y _|_ to peptide plane  ="
print "=                                                          ="
print "= 3) Z || to C(O)-N bond, Y _|_ to peptide plane           ="
print "=                                                          ="
print "=                                                          ="
print "=                    C(a)'    H                            ="
print "=                     \      /                             ="
print "=                      C----N                              ="
print "=                     //     \                             ="
print "=                    O        C(a)                         ="
print "=                                                          ="
print "============================================================"
print ""
#######################
# CHECK INPUT OPTIONS #
#######################
argc = len(sys.argv)
if (argc < 4):
	print ">>> ERROR <<<"
	print "Some arguments are missing. Please, run as: python NH_tilt.py f1 f2 method"
	print "where:"
	print "* f1 : xyz ball-and-sticks file of protein, generated with Babel 1.1 (-obs output option)"
	print "* f2 : xyz file generated by DiTe, containing protein coordinates expressed in diffusion frame (MF)"
	print "* method : integer [1, 2 or 3] specifying how to calculate the director frame."
	print ""
	exit(1)
connectivityFile = sys.argv[1]
coordinatesFile = sys.argv[2]
directorFrame = int(sys.argv[3])
#############################
# SOME INITIAL DECLARATIONS #
#############################
# list of atomic types
atomType=[]
# list of XYZ coordinates
xyz=[[],[],[]]
# list of connectivity
connect=[]
# file string
fstring=""
######################
# LEARN CONNECTIVITY #
######################
f = open(connectivityFile,"r")
# first two lines
fstring = f.readline()
fstring = f.readline()
# start reading file
fstring = f.readline()
while (fstring!="" and fstring!="\n"):
	splittedString = fstring.split()
	# store atom type
        atomType.append(splittedString[0])
	# store connectivity
	connect.append(splittedString[4:len(splittedString)])
	# read next line
	fstring = f.readline()
f.close()
###################
# NUMBER OF ATOMS #
###################
natoms = len(atomType)
####################################
# READ COORDINATES EXPRESSED IN MF #
####################################
f = open(coordinatesFile,"r")
# first line
fstring = f.readline()
# start reading geometry
fstring = f.readline()
while (fstring!="" and fstring!="\n"):
	splittedString = fstring.split()
	# store Cartesian coordinates
	xyz[0].append(float(splittedString[1]))
	xyz[1].append(float(splittedString[2]))
	xyz[2].append(float(splittedString[3]))
	# read next line
	fstring = f.readline()
f.close()
#**********************************************************
# DEBUGGING INFO : echo geometry file in MF
#**********************************************************
#for i in range(0,len(atomType)):
#	print atomType[i], xyz[0][i], xyz[1][i], xyz[2][i],
#	for j in range(0,len(connect[i])):
#		print connect[i][j],
#	print ""
#*********************************************************
####################################################
# DETERMINE EULER ANGLES MF -> VF FOR EACH RESIDUE #
####################################################
# number of residues
nres = 1
# list of residues
res=[]
# N coord ----------------- N
ncoord = [0.0,0.0,0.0]
# C0 coord ---------------- C(a)'
c0coord = [0.0,0.0,0.0]
# C1 coord ---------------- C
c1coord = [0.0,0.0,0.0]
# C2 coord ---------------- C(a)
c2coord = [0.0,0.0,0.0]
# lists of director cosines matrices
rV = []
rD = []
# lists of origin of VF frames
OX = []
OY = []
OZ = []
# dictionary of functions to use based on method
VF = {1:Z_NH, 2:Z_CC, 3: Z_CN}#, 4: OmegaD}
#                                        #
# SEARCH FOR C0-C1(=O)-N(H)-C2 STRUCTURE #
#                                        # 
for ia1, a1 in enumerate(atomType):
	if (a1=='N' or a1=='n'):
		# search for Prolin residue: no N-H, but counts as a residue
		if (len(connect[ia1]) == 3):
			ia2 = int(connect[ia1][0])-1
			ia3 = int(connect[ia1][1])-1
			ia4 = int(connect[ia1][2])-1
			atomsList = atomType[ia2]+atomType[ia3]+atomType[ia4]
			if (atomsList == "CCC" or atomsList == "ccc"):
				nres = nres + 1
		# search for other residues
		if (len(connect[ia1]) == 2):
			ia2 = int(connect[ia1][0])-1
			ia3 = int(connect[ia1][1])-1
			atomsList = atomType[ia2] + atomType[ia3]
			Nfound = 0
			# check if atoms bonded to N are both C atoms 
			if (atomsList == "CC" or atomsList == "cc" or atomsList == "Cc" or atomsList == "cC"):
				# find if and wich of the C atoms is carboxylic
				for a4 in connect[ia2]:
					ia4 = int(a4) - 1
					if (atomType[ia4] == "O" or atomType[ia4] == "o"): 
						Nfound = 1
						nres = nres + 1
				for a5 in connect[ia3]:
					ia5 = int(a5) - 1
					if (atomType[ia5] == "O" or atomType[ia5] == "o"):
						Nfound = 1
						nres = nres + 1
						ia2, ia3 = ia3, ia2
				# locate the C(i-1) bonded to C1
				for a6 in connect[ia2]:
					ia6 = int(a6) - 1
					if (atomType[ia6] == "C" or atomType[ia6] == "c"):
						ia4 = ia6
						break
			#                                                                      #
			# CALCULATE THE DIRECTOR COSINES IF THE BACKBONE N ATOM HAS BEEN FOUND #
			#                                                                      #
			if (Nfound > 0):
				# update residues array
				res.append(nres)
				#**********************************************************
				# DEBUGGING INFO : residue number and C1-N-C2 atoms numbers
				#**********************************************************
				# print nres
				# print ia1+1, ncoord
				# print ia2+1, c1coord
				# print ia3+1, c2coord
				#**********************************************************
				# assign atoms coordinates to vectors ncoord, c1coord and c2coord
				ncoord[0], ncoord[1], ncoord[2] = xyz[0][ia1], xyz[1][ia1], xyz[2][ia1]
				c1coord[0], c1coord[1], c1coord[2] = xyz[0][ia2], xyz[1][ia2], xyz[2][ia2]
				c2coord[0], c2coord[1], c2coord[2] = xyz[0][ia3], xyz[1][ia3], xyz[2][ia3]
				c0coord[0], c0coord[1], c0coord[2] = xyz[0][ia4], xyz[1][ia4], xyz[2][ia4]
				# Calculate ROTATION matrices
				R, x0, y0, z0 = VF.get(directorFrame,UNKNOWN_METHOD)(ncoord,c0coord,c1coord,c2coord)
				rV.append(R)
				R, xD, yD, zD = VF.get(1,UNKNOWN_METHOD)(ncoord,c0coord,c1coord,c2coord)
				rD.append(R)
				# append VF origin
				OX.append(x0)
				OY.append(y0)
				OZ.append(z0)
# reset nres to the real number of residues
nres = len(res);
##########################
# CALCULATE EULER ANGLES #
##########################
tmplist = coordinatesFile.split(".xyz")
eulerFile = tmplist[0] + ".eul"
f1 = open(eulerFile,"w")
eulV=[[],[],[]]
eulD=[[],[],[]]
X, Y, Z = 0, 1, 2
rad2deg = 180.0 / pi
ZERO = 1.0e-10
for i in range (0,len(res)):
	RPV = rV[i]
	RPD = rD[i]
	RVD = numpy.dot(RPD,RPV.T)
	# Euler matrix is the transpose of the rotation matrix
	RPV = RPV.T
	RVD = RVD.T
	for m in range(0,3):
		for n in range (0,3):
			if (abs(RPV[m,n] - 0.0) < ZERO) : RPV[m,n] =  0.0
			if (abs(RVD[m,n] - 0.0) < ZERO) : RVD[m,n] =  0.0
			if (abs(RPV[m,n] - 1.0) < ZERO) : RPV[m,n] =  1.0
			if (abs(RVD[m,n] - 1.0) < ZERO) : RVD[m,n] =  1.0
			if (abs(RPV[m,n] + 1.0) < ZERO) : RPV[m,n] = -1.0
			if (abs(RVD[m,n] + 1.0) < ZERO) : RVD[m,n] = -1.0
	# Euler angles from VF to MF
	RPV = RPV.T
	eulV[0].append(atan2(RPV[Y,Z],RPV[X,Z])*rad2deg)
	eulV[1].append(atan2(sqrt(1.0-RPV[Z,Z]*RPV[Z,Z]),RPV[Z,Z])*rad2deg)
	eulV[2].append(atan2(RPV[Z,Y],-RPV[Z,X])*rad2deg)
	# Euler angles from OF to DF
	eulD[0].append(atan2(RVD[Y,Z],RVD[X,Z])*rad2deg)
	eulD[1].append(atan2(sqrt(1.0-RVD[Z,Z]*RVD[Z,Z]),RVD[Z,Z])*rad2deg)
	eulD[2].append(atan2(RVD[Z,Y],-RVD[Z,X])*rad2deg)
	f1.write("%3d\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f\t%3.3f"%(res[i],eulV[0][i],eulV[1][i],eulV[2][i],eulD[0][i],eulD[1][i],eulD[2][i]) + "\n")
f1.close()
####################
# OUTPUT FOR PYMOL #
####################
red   = ",0.1,1.0,0.0,0.0,1.0,0.0,0.0]"
green = ",0.1,0.0,1.0,0.0,0.0,1.0,0.0]"
blue  = ",0.1,0.0,0.0,1.0,0.0,0.0,1.0]"
tmplist = coordinatesFile.split(".xyz")
pymolFile = tmplist[0] + ".py"
f1 = open(pymolFile,"w")
f1.write("from pymol import *" + "\n")
f1.write("from pymol.cgo import *" + "\n")
f1.write("cmd.load ('"+coordinatesFile + "')" + "\n")
for i in range (0,len(OX)):
	R = rV[i]
	str0 = "=[CYLINDER," + str(OX[i]) + "," + str(OY[i]) + "," + str(OZ[i]) + ","
	xx, yy, zz = OX[i]+R[X,X], OY[i]+R[X,Y], OZ[i]+R[X,Z]
	f1.write("x" + str(i) + str0 + str(xx) + "," + str(yy) + "," + str(zz) + red + "\n")
	f1.write("cmd.load_cgo(x" + str(i) + ",'x" + str(i) + "')" + "\n")
	xx, yy, zz = OX[i]+R[Y,X], OY[i]+R[Y,Y], OZ[i]+R[Y,Z]
	f1.write("y" + str(i) + str0 + str(xx) + "," + str(yy) + "," + str(zz) + green + "\n")
	f1.write("cmd.load_cgo(y" + str(i) + ",'y" + str(i) + "')" + "\n")
	xx, yy, zz = OX[i]+R[Z,X], OY[i]+R[Z,Y], OZ[i]+R[Z,Z]
	f1.write("z" + str(i) + str0 + str(xx) + "," + str(yy) + "," + str(zz) + blue + "\n")
	f1.write("cmd.load_cgo(z" + str(i) + ",'z" + str(i) + "')" + "\n")
f1.close()

