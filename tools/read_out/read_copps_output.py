import sys, os
###############
# CHECK INPUT #
###############
argc = len(sys.argv)
inFileName = ""
outFileName = "fit_results.dat"
if (argc < 2):
	print ""
	print "Usage: python read_copps_output.py fin fout"
	print "where:"
	print "* fin [required] is the input file name (usually XXX_copps.output)"
	print "* fout [optional] is the output file name (default: fit_results.dat)"
	print ""
	exit(1)
else:
	inFileName = sys.argv[1]
	if (argc > 2):
		outFileName = sys.argv[2]
#####################
# PREPARE THE FILES #
#####################
fin = open(inFileName,"r")
fout = open(outFileName,"w")
########################
# SOME INITIALIZATIONS #
########################
nComp = 0
srlsData = []
nmrData = []
nRelaxData = 0
relaxData = []
nFields = 0
fields=[]
nFitData = 0
fitData = []
fitDataUnits = []
tableTitles = []
chisq = 0.0
nDoF = 0
##############################
# LOOP OVER INPUT FILE LINES #
##############################
fileLine = fin.readline()
while (len(fileLine) > 0):
	## 1. Locate the '* Parameters' line
	if (fileLine == "* Parameters:\n"):
		fileLine = fin.readline() # empty line
		fileLine = fin.readline()
		## 2. Store SRSL fitting parameters
		while (fileLine != "\n"):
			srlsData.append(fileLine.replace("\n",""))
			fileLine = fin.readline()
		fileLine = fin.readline() # * Relaxation data
		fileLine = fin.readline() # empty line
		fileLine = fin.readline() # table titles
		fileLine = fin.readline() # table line
		fileLine = fin.readline() # empty line
		fileLine = fin.readline()
		## 3. Store NMR data
		while (fileLine[0] != "*"):
			splitLine = fileLine.split()
			if (len(splitLine) == 5):
				nmrData.append(fileLine.replace("\n",""))
			fileLine = fin.readline()
		## 4. Store Chi square
		splitLine = fileLine.split()
		chisq = float(splitLine[4])
		## 5. Order parameters
		fileLine = fin.readline() # empty line
		fileLine = fin.readline() # Order parameters
		fileLine = fin.readline() # empty line
		fileLine = fin.readline() # S20
		srlsData.append(fileLine.replace("\n",""))
		fileLine = fin.readline() # S22
		srlsData.append(fileLine.replace("\n",""))
		## 6. Determine if this was the last fitting step for component nComp
		fileLine = fin.readline() # Sxx
		fileLine = fin.readline() # Syy
		fileLine = fin.readline() # Szz
		fileLine = fin.readline() # empty line
		fileLine = fin.readline() # empty line
		fileLine = fin.readline() # * Elapsed time
		fileLine = fin.readline() # empty line
		fileLine = fin.readline()
		#print fileLine
		splitLine = fileLine.split()
		if (len(splitLine) == 6):
			## 7. If this is the first component, print table titles
			if (nComp == 0):
				### 7.A SRLS
				nFitData = len(srlsData) - 2
				for i in range (0,nFitData):
					splitLine = srlsData[i].split("=")
					fitData.append(splitLine[0])
					splitLine = splitLine[1].split()
					fitDataUnits.append(splitLine[1])
				srlsLine = [" " for i in range (0,len(srlsData)+1)]
				splitLine = nmrData[0].split()
				### 7.B NMR
				refStr = splitLine[0]
				relaxData.append(refStr[0:refStr.index(":")])
				for i in range (1,len(nmrData)):
					splitLine = nmrData[i].split()
					if (splitLine[0] == refStr):
						break
					tmpStr = splitLine[0]
					relaxData.append(tmpStr[0:tmpStr.index(":")])
				nRelaxData = len(relaxData)
				nFields = len(nmrData) / nRelaxData
				for i in range (0, len(nmrData), nRelaxData):
					splitLine = nmrData[i].split()
					fields.append(splitLine[1])
				nDoF = 1 # Because C++OPPS v >= 1.1.1 just gives the reduced chi-square
				### 7.C write output file heading
				fout.write("# Input file : " + inFileName + "\n")
				fout.write("# Number of fitting parameters : " + str(nFitData) + "\n")
				fout.write("# Fitting parameters :" + "\n")
				for i in range (0, nFitData):
					fout.write("#  * " + fitData[i] + " / " + fitDataUnits[i] + "\n")
				fout.write("# Number of frequencies : " + str(nFields) + "\n")
				fout.write("# Frequencies : " + "\n")
				for i in fields:
					fout.write("# * " + i + " MHz" + "\n")
				fout.write("# Number of NMR data per frequency : " + str(nRelaxData) + "\n")
				fout.write("# NMR data : " + "\n")
				for i in relaxData:
					fout.write("# * " + i + "\n")
				fout.write("#\n")
				fout.write("# Table will show data in the following order:" + "\n")
				fout.write("# Component | SRLS data | S20 | S22 | Red. Chi square | NMR data freq 1 | ... | NMR data field nFreqs" + "\n")
				fout.write("#\n")
				fout.write("# Notes : " + "\n")
				fout.write("# * Red. Chi square is Chi square divided by nDegreesOfFreedom = nRelaxData - nFitData - 1" + "\n")
				fout.write("# * NMR data represented as: theoretical, experimental, %error" + "\n")
				fout.write("#\n")
				fout.write("#\n")
			## 8. Print data to fout
			fout.write(str(nComp+1) + "\t")
			for i in srlsData:
				spliti = i.split("=")
				spliti = spliti[1].split()
				fout.write(spliti[0] + "\t")
			fout.write(str(chisq/float(nDoF)) + "\t")
			for i in nmrData:
				spliti = i.split()
				fout.write(spliti[2] + "\t")
				fout.write(spliti[3] + "\t")
				fout.write(spliti[4] + "\t")
			fout.write("\n")
			## 9. Update components counter
			nComp = nComp + 1
	## 10. reset data
	del nmrData[:]
	del srlsData[:]
	fileLine = fin.readline()
fin.close()
fout.close()
