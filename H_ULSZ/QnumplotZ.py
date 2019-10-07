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
from scipy.interpolate import interp1d
pi = 3.14159265



##### ----------------------
def main():
	
	filename="Sz_ULSZ.dat"
	file = open(filename,'r')
	lines = file.readlines()
	file.close()
	
	filename="QnumZ.sh"
	file = open(filename,'r')
	lines2 = file.readlines()
	file.close()

	lines_lam = lines2[2]
	lines_lam = lines_lam.replace("for Hz in ","")
	
	Nl = len(lines_lam.split(" "))
	Sz = np.zeros(Nl)
	Lam = np.zeros(Nl)

	for i in range(0,Nl):
		temp = lines_lam.split(" ")	
		out = float(temp[i])
		Lam[i] = out

	for i in range(0,Nl):
		temp = lines[0].split(" ")
		out = float(temp[i])
		Sz[i] = out	
	
	for i in range(0,Nl):
		temp = lines_lam.split(" ")	
		out = float(temp[i])
		Lam[i] = out
	#print(Lam)


	
	
	#----------------------------plot--------------------------

	plt.rc('font', family='serif')	
	fig = plt.figure(figsize=(10,4)) 
	ax=fig.add_subplot(111)

	markersize=5
	plt.xlabel('h',fontsize=24)
	plt.xticks(fontsize=18)
	plt.yticks(fontsize=18)
	plt.plot(Lam, Sz, 'o-',linewidth=2,color='black',markersize=markersize,label = 'ULSZ')

	plt.tick_params(
	    axis=u'y',          # changes apply to the y-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='on',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='on', right='off', left='on', labelleft='on') # labels along the bottom edge are off

	plt.tick_params(
	    axis=u'x',          # changes apply to the x-axis
	    which='both',      # both major and minor ticks are affected
	    bottom='on',      # ticks along the bottom edge are off
	    top='off',         # ticks along the top edge are off
	    labelbottom='on', right='off', left='on', labelleft='on') # labels along the bottom edge are off

	plt.grid(which='major',linestyle='--',alpha=0.25)
	plt.tick_params(which=u'major',axis=u'both',length=6,color="black",direction='in')
	plt.tick_params(which=u'minor',axis=u'x',length=3,color="black",direction='in')
	plt.tick_params(which=u'minor',axis=u'y',length=3,color="black",direction='in')
	plt.locator_params(nbins=4)
	minor_locator = tkr.AutoMinorLocator(4)
	ax.xaxis.set_minor_locator(minor_locator)
	minor_locator = tkr.AutoMinorLocator(2)
	ax.yaxis.set_minor_locator(minor_locator)
	plt.legend(loc='best',ncol=1, fontsize=18, frameon=False,columnspacing=1,handlelength=1.0,handletextpad=0.4)
	plt.savefig('Sz_ULSZ.pdf',bbox_inches='tight',dpi=300)
	plt.close()


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
	orb=2
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
			BB[i,j] = mat[ib,jb]
			AB[i,j] = mat[ia,jb]

	return AA, BB, AB

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
			if(i==j): 
				m[i,j] = 2.0;
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



def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
##### ----------------------



if __name__ == '__main__':
	sys.argv ## get the input argument
	total = len(sys.argv)
	cmdargs = sys.argv	
	main()
































