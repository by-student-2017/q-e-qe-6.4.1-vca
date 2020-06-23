#!/usr/bin/python
import sys,math

# pwout2in.py
# made by mkanzaki@me.com
# Utility code for Quantum-Espsresso.
# This code converts pw.x output (optimization) to new input file.
# Usage:
# >./pwout2in.py test
# Then test.new.in will be produced.
# test.in and test.out must exist.
#

cell = ['' for i in range(4)]
position = ['' for i in range(1000)]

# conversion for radian 
rd = math.pi/180.0
# conversion factor for atomic unit to angstrom 
bohr = 0.5291772108

# reading output file
# check xtl file name is given or not
if len(sys.argv) == 1:
	print 'No file name provided!'
	print 'for example: > ./pwout2in.py test'
	print 'test.out and test.in must exist.'
	exit()
try:
	f1 = sys.argv[1] + '.in'
	file1 = open(f1,'r')
	f2 = sys.argv[1] + '.out'
	file2 = open(f2,'r')
except IOError, (errno, msg):
	print 'File open error!'
	exit()
# Get file name
f3 = sys.argv[1] + '.new.in'
#f3 = 'temp.in'
# Open new file .in
file3 = open(f3,'w')

# find a line containing "Final" in .out file
while True:
	ftext = file2.readline()
        if 'Final' in ftext :
                break
	if ftext=="":
		print 'Final optimization result was not found in this output file!'
		file1.close()
		file2.close()
		file3.close()
		exit()		
# find cell parameters
while True:
	ftext = file2.readline()
	if 'CELL_PARAMETERS' in ftext :
		tmp = ftext.strip()
		break
# get alat 
s = tmp.index('alat=')+5
t = tmp.index(')')
alat = float(tmp[s:t])
# reading cell matrix
cell[1] = file2.readline()
cell[2] = file2.readline()
cell[3] = file2.readline()
# find atomic positions in .out file
while True:
	ftext = file2.readline()
	if 'ATOMIC_POSITIONS' in ftext :
		break
# read atomic positions in .out file
i = 1
while True:
	ftext = file2.readline()
	if 'End final coordinates' in ftext :
		break
	else:
		position[i] = ftext
		i = i + 1
		#print str(i)
natom = i
file2.close()

# read original .in file
while True:
	ftext = file1.readline()
	if 'ibrav' in ftext :
		break
	file3.write(ftext)
file3.write('    ibrav = 0 ,\n')
file3.write('    celldm(1) = ' + str(alat) + ' ,\n')
while True: # ignore celldm(i) i>1
	ftext = file1.readline()
	if 'nat' in ftext :
		break
file3.write(ftext)
while True:
        ftext = file1.readline()
        if "&cell" in ftext :
                file3.write(ftext)
                break
	file3.write(ftext)
while True:
	ftext = file1.readline()
	if '/' in ftext :
		file3.write(ftext)
		break
	file3.write(ftext)
#
ftext = file1.readline()
if 'CELL_PARAMETERS' in ftext :
		tmp = file1.readline()
		tmp = file1.readline()
		tmp = file1.readline()
		tmp = file1.readline()
# write new cell matrix
file3.write('CELL_PARAMETERS (alat= ' + str(alat) + ')\n')
file3.write(cell[1])
file3.write(cell[2])
file3.write(cell[3])
file3.write('ATOMIC_SPECIES\n')		
while True:
	ftext = file1.readline()
	if 'ATOMIC_POSITIONS' in ftext :
		break
	file3.write(ftext)
file3.write(ftext)
for i in range(1,natom): # write new positions
	file3.write(position[i])
while True:
	ftext = file1.readline()
	if 'K_POINTS' in ftext :
		break
file3.write(ftext)
ftext = file1.readline()
file3.write(ftext)
file1.close()
file3.close()

