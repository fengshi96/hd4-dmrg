import sys, re, math, random		## for passing an argument and list of variables ## regexes, math functions, random numbers
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.interpolate
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.ticker as tkr
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from matplotlib.patches import Arc
from subprocess import call
pi = 3.14159265



##### ----------------------
def main(total, cmdargs):
	if total != 3:
			print (" ".join(str(x) for x in cmdargs))
			raise ValueError('I did not ask for arguments')
	Ns=int(cmdargs[1])
	Lambda=str(cmdargs[2])
	Orbitals=1
	NumberOfOrb=Ns*Orbitals;
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	matplotlib.rcParams['axes.linewidth'] = 2 #set the value globally


	### -------- onsite observables -------------------
	filename="Observables_Local2.dat"	
	string="j Lx Ly Lz = "
	Lx = np.zeros(Ns); Ly = np.zeros(Ns); Lz = np.zeros(Ns); 


	file = open(filename,'r')
	lines = file.readlines()
	file.close()
	for i in range(0,len(lines)):
		line = lines[i]
		
		if re.search(string, line, re.I):
			matchindex=i+1
			break
	try:
		matchindex
	except NameError:
		print ("***match string", string)
		print ('***MYERROR:matchindex in read_density function is not defined - no match found')

	for i in range(0,Ns):
		idx = matchindex+i
		#print()	
		line = lines[idx];
		temp = line.split(" ")
		out0 = ConvertImag(temp[1]); out1 = ConvertImag(temp[2]); out2 = ConvertImag(temp[3]); 
		Lx[i] = out0.real; Ly[i] = out1.real; Lz[i] = out2.real;


	Lz_avg = Lz.sum()/Ns
	
	outfileS = open("../../Sz_ULSZ.dat","a")
	outfileS.write(str(Lz_avg)+" ")
	outfileS.close()














#### ========================================================================================
#### ========================= Functions ====================================================
#### ========================================================================================
#### ========================================================================================
#### ========================= Functions ====================================================
#### ========================================================================================
def read_density(NumberOfSites,str,matchstring):
	file = open(str,'r')
	lines = file.readlines()
	file.close()
	for i in range(0,len(lines)):
			line = lines[i]
			if re.search(matchstring, line, re.I):
				matchindex=i+1
				break
	try:
		matchindex
	except NameError:
		print ("***match string", matchstring)
		print ('***MYERROR:matchindex in read_density function is not defined - no match found')
	rows = NumberOfSites
	cols = 2
	m = [[0.0 for x in range(cols)] for y in range(rows)] ## allocate m list
	for i in range(0,rows):
		line = lines[matchindex+1+i];
		temp = line.split(" ")
		val = ConvertImag(temp[1])
		m[i][0] = val.real;
		m[i][1] = val.imag;  	### defined 2D Matrix
	return np.asarray(m);

##### ----------------------
def split(mat):
	orb=1
	nrows = mat.shape[0];
	ncols = mat.shape[1]; ##if mat else 0
	mat = np.asarray(mat)
	assert(nrows==ncols)
	fatrows=int(nrows/orb);
	fatcols=int(ncols/orb);

	AA = np.zeros((fatrows,fatcols));
	BB = np.zeros((fatrows,fatcols));
	AB = np.zeros((fatrows,fatcols));

	for i in range(0,fatrows):
		ia = i*orb+0
		ib = i*orb+1
		for j in range(0,fatcols):
			ja = j*orb+0
			jb = j*orb+1
	
			AA[i,j] = mat[ia,ja]


	return AA

##### ----------------------
def ConvertImag(s):
	repart = float(s.split(",")[0].split("(")[1])
	impart = float(s.split(",")[1].split(")")[0])
	return complex(repart,impart)

##### ----------------------
def RepresentsInt(s):
    try:
        int(s)
        return True
    except ValueError:
        return False

##### ----------------------
def readMatrix(str):
	file = open(str,'r')
	lines = file.readlines()
	file.close()

	Type = lines[0].split(" ")[2]
	Type = ConvertImag(Type);
	rows = int(Type.real)
	cols = int(Type.imag)

	counter = 0
	m = np.zeros((rows,cols)) ### [[0.0 for x in range(rows)] for y in range (cols)]
	for i in range(1,rows+1):
		line = lines[i];
		line = " ".join(line.split())
		#print(line)
		temp = line.split(" ")
		for j in range(0,cols):
			#m[j+i*cols] = float(temp[j])  ### defined 1D aray
			val = ConvertImag(temp[j]);
			m[i-1,j] = val.real	   ### defined 2D Matrix
		counter = counter + 1
	file.close()

	for i in range(0,rows):
		for j in range(i,cols):
			m[j,i] = m[i,j]	   ### defined 2D Matrix
	return m; 
##### ----------------------

	
def PrintMatrix(m):
	string="""
DegreesOfFreedom=1
GeometryKind=LongRange
GeometryOptions=none
Connectors"""
	
	num_rows = len(m);
	num_cols = m.shape[1] #len(m[0]) if m else 0
	string += " "+str(num_rows)+" "+str(num_cols)
	print (string )
	for i in range(0,num_rows):
		print (" ".join(str(x) for x in m[i]))
##### ----------------------


if __name__ == '__main__':
	sys.argv ## get the input argument
	total = len(sys.argv)
	cmdargs = sys.argv	
	main(total, cmdargs)
































