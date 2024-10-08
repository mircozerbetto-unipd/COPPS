import sys
import string
from string import *
#######################
##     HELLO MSG     ##
#######################
print ""
print "================================"
print "TXT TO EXP CONVERTER FOR C++OPPS"
print "================================"
print ""
#######################
##    INPUT CHECK    ##
#######################
argc = len(sys.argv);
print "Important options:"
print "  -i <inputFileName>"
print "  -o <outputFileName>"
print ""
print "Other options:"
print "  -r1r2 R/T [default R]"
print "  -u Hz/s/ms [default Hz for R, s for T]"
print "  -f field/freq [default field]"
print ""
if (argc < 5):
	print ">>> ERROR <<<"
	print "Some important options are missing"
	print ""
	exit(1)
#######################
##  ASSIGN OPTIONS   ##
#######################
# check input/output file names
if (sys.argv[1] == "-i"):
	inputFileName = sys.argv[2];
	if (sys.argv[3] == "-o"):
		outputFileName = sys.argv[4]
	elif (sys.argv[3] == "-i"):
		print ">>> ERROR <<<"
		print "Duplicate -i option"
		print ""
		exit(1)
	else:
		print ">>> ERROR <<<"
		print "Unrecognized option "+sys.argv[3];
		print ""
		exit(1)
elif (sys.argv[1] == "-o"):
	outputFileName = sys.argv[2];
	if (sys.argv[3] == "-i"):
		inputFileName = sys.argv[4]
	elif (sys.argv[3] == "-o"):
		print ">>> ERROR <<<"
		print "Duplicate -o option"
		print ""
		exit(1)
	else:
		print ">>> ERROR <<<"
		print "Unrecognized option "+sys.argv[3];
		print ""
		exit(1)
else:
	print ">>> ERROR <<<"
	print "Unrecognized option "+sys.argv[1];
	print ""
	exit(1)
# set some defaults
r1r2 = "r"
r1r2Found = False
units = "hz"
unitsFound = False
field = "field"
fieldFound = False
# override defaults based on user choices
for i in range (3,argc):
	if (sys.argv[i] == "-r1r2"):
		if (r1r2Found == True):
			print ">>> WARNING <<<"
			print "duplicated option -r1r2 found on argument " + str(i)
		r1r2Found = True
		i = i+1
		r1r2 = lower(sys.argv[i])
		if (r1r2 == "t" and unitsFound == False):
			units = "s"
		if (r1r2 != "t" and r1r2 != "r"):
			print ">>> ERROR <<<"
			print "Argument " + r1r2 + " not recognized for oprion -r1r2"
			print ""
			exit(1)
	elif (sys.argv[i] == "-u"):
		if (unitsFound == True):
			print ">>> WARNING <<<"
			print "duplicated option -u found on argument " + str(i)
		unitsFound = True
		i = i+1
		units = lower(sys.argv[i])
		if (units != "hz" and units != "s" and units != "ms"):
			print ">>> ERROR <<<"
			print "Argument " + units + " not recognized for oprion -u"
			print ""
			exit(1)
	elif (sys.argv[i] == "-f"):
		if (fieldFound == True):
			print ">>> WARNING <<<"
			print "duplicated option -f found on argument" + str(i)
		fieldFound = True
		i = i+1
		field = lower(sys.argv[i])
		if (field != "freq" and field != "field"):
			print ">>> ERROR <<<"
			print "Argument " + field + " not recognized for oprion -f"
			print ""
			exit(1)
# check units
if (r1r2 == "r" and units != "hz"):
	print ">>> ERROR <<<"
	print "Incompatible units for relaxation rates (R1 and R2). The '-r1r2 R' option is compatible only with the '-u Hz' option"
	print ""
	exit(1)
if (r1r2 == "t" and units == "hz"):
	print ">>> ERROR <<<"
	print "Incompatible units for relaxation times (T1 and T2). The '-r1r2 T' option is compatible only with the '-u s' or '-u ms' options"
	print ""
	exit(1)
#######################
##     CONSTANTS     ##
#######################
gammaH = 42.576 # MHz/T
#######################
##    READ INPUT     ##
#######################
f1 = open(inputFileName);
flines = f1.readlines();
f1.close()
# echo infput file
print "Echoing input file content:"
print ""
print "--------------- "+inputFileName+" --------------- START"
for line in flines:
	print line
print "--------------- "+inputFileName+" --------------- EOF"
print ""
# number of fields / frequencies
fields = flines[0].split()
nFields = len(fields)
print "================================"
print ""
print str(nFields) +" fields found in input:"
for i in range (0,nFields):
	if (field == "field"):
		fields[i] = float(fields[i])
	else:
		fields[i] = float(fields[i])/gammaH
	print "  " + str(fields[i]) + " T -> " + str(float(fields[i])*gammaH) + " MHz"
# read data
nData = 1 + 6*nFields
listOfResidues = []
nComponents = 0
f1 = open(outputFileName,"w")
for i in range (1,len(flines)):
	dataLine = flines[i].split()
	if (len(dataLine) > 1):
		if (len(dataLine) != nData):
			print ">>> ERROR <<<"
			print "Wrong lenght of row " + str(i+1) + " in input file. Found " + str(len(dataLine)) + " numbers instead of " + str(nData)
			print ""
			exit(1)
		nComponents = nComponents + 1
		listOfResidues.append(dataLine[0])
		f1.write( "# data for component " + str(nComponents) + " - RESIDUE " + str(dataLine[0]) + "\n" )
		f1.write( "Component:" + str(nComponents) + "\n" )
		f1.write( "Residue:" + str(dataLine[0]) + "\n" )
		for j in range (0,nFields):
			step = j*6
			if (r1r2 == "r"):
				units = "s"
				dataLine[1+step] = str( 1.0/float(dataLine[1+step]) )
				dataLine[2+step] = str( float(dataLine[2+step]) * float(dataLine[1+step]) * float(dataLine[1+step]) )
				dataLine[3+step] = str(1.0/float(dataLine[3+step]))
				dataLine[4+step] = str( float(dataLine[4+step]) * float(dataLine[3+step]) * float(dataLine[3+step]) )
			f1.write( "# filed " + str(j+1) + " - " + str(fields[j]) + " T" + "\n" )
			f1.write( "Field:"+str(fields[j]*gammaH)+":MHz" + "\n")
			f1.write( "T1:"+str(round(float(dataLine[1+step]),4))+":"+str(round(float(dataLine[2+step]),4))+":"+units + "\n" )
			f1.write( "T2:"+str(round(float(dataLine[3+step]),4))+":"+str(round(float(dataLine[4+step]),4))+":"+units + "\n" )
			f1.write( "NOE:"+str(round(float(dataLine[5+step]),4))+":"+str(round(float(dataLine[6+step]),4))+":-" + "\n" )
		f1.write( "# M1F -> VF Euler angles" + "\n" )
		f1.write( "Euler_angles:0.0:0.0:0.0:deg" + "\n" )
f1.close()
print ""
print "================================"
print ""
print str(nComponents) + " components found. Echoing Component-Residue map:"
print ""
for i in range (0,nComponents):
	print str(i+1) + "\t" + str(listOfResidues[i])
