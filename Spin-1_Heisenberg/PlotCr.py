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
	Orbitals=2
	NumberOfOrb=Ns*Orbitals;
	plt.rc('text', usetex=True)
	plt.rc('font', family='serif')
	matplotlib.rcParams['axes.linewidth'] = 2 #set the value globally


	### -------- onsite observables -------------------
	filename="Observables_Local.dat"	
	string="j Sx Sy Sz Lx Ly Lz "
	Sx = np.zeros(Ns); Sy = np.zeros(Ns); Sz = np.zeros(Ns); 
	Lx = np.zeros(Ns); Ly = np.zeros(Ns); Lz = np.zeros(Ns); 
	Jz = np.zeros(Ns); 


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
		line = lines[idx];
		temp = line.split(" ")
		out0 = ConvertImag(temp[1]); out1 = ConvertImag(temp[2]); out2 = ConvertImag(temp[3]); 
		out3 = ConvertImag(temp[4]); out4 = ConvertImag(temp[5]); out5 = ConvertImag(temp[6]); 
		out6 = ConvertImag(temp[7]); 

		Sx[i] = out0.real; Sy[i] = out1.real; Sz[i] = out2.real;
		Lx[i] = out3.real; Ly[i] = out4.real; Lz[i] = out5.real;
		Jz[i] = out6.real;





	### -------- two-point correlations -------------------
	filename="Observables_SxSx.dat"	
	CxCx = readMatrix(filename);
	filename="Observables_SySy.dat"	
	CyCy = readMatrix(filename);
	filename="Observables_SzSz.dat"	
	CzCz = readMatrix(filename);

	Cx_ss,Cx_ll,Cx_sl = split(CxCx);
	Cy_ss,Cy_ll,Cy_sl = split(CyCy);
	Cz_ss,Cz_ll,Cz_sl = split(CzCz);
	StSt = Cx_ss + Cy_ss + Cz_ss;
	LtLt = Cx_ll + Cy_ll + Cz_ll;
	StLt = Cx_sl + Cy_sl + Cz_sl;


	for i in range(0,Ns):
		StSt[i,i] = 2.0; 
		LtLt[i,i] = 2.0; 
	JtJt=StSt+LtLt+StLt+StLt;

	#print StSt

	for i in range(0, Ns):
		for j in range(0, Ns):
			testss = Sx[i] * Sx[j] + Sy[i] * Sy[j] + Sz[i] * Sz[j]
			StSt[i,j] = StSt[i,j] - testss

			testll = Lx[i] * Lx[j] + Ly[i] * Ly[j] + Lz[i] * Lz[j]
			LtLt[i,j] = LtLt[i,j] - testll

			testsl = Sx[i] * Lx[j] + Sy[i] * Ly[j] + Sz[i] * Lz[j]
			StLt[i,j] = StLt[i,j] - testsl

			testsl_T = Sx[j] * Lx[i] + Sy[j] * Ly[i] + Sz[j] * Lz[i]
			test = testss + testll + testsl + testsl_T
			JtJt[i,j] = JtJt[i,j] - test

	#print StSt; 

	Sk = np.zeros(Ns+1);
	Lk = np.zeros(Ns+1);
	SLk = np.zeros(Ns+1);
	Jk = np.zeros(Ns+1);
	for k in range(0,Ns+1):
		kx = float(k)*2.0*pi/Ns
		for i in range(0,Ns):
			for j in range(0,Ns):
				Sk[k] += StSt[i,j]*np.cos(kx*(i-j))
				Lk[k] += LtLt[i,j]*np.cos(kx*(i-j))
				SLk[k] += StLt[i,j]*np.cos(kx*(i-j))
###################
				Jk[k] += JtJt[i,j]*np.cos(kx*(i-j))


	# -- Fourier transform to momentum space -- 
	filename4="SkNk.dat"
	targetfile4 = open(filename4, 'w')
	for k in range(0,Ns+1):
		kx = float(k)*2.0*3.14159265/Ns
		string = str(kx/pi) + " " + str(Sk[k]/(Ns*Ns))
		string += " " + str(Lk[k]/(Ns*Ns))
		string += " " + str(SLk[k]/(Ns*Ns))
		string += " " + str(Jk[k]/(Ns*Ns))
		#print string
		targetfile4.write(string)
		targetfile4.write(' \n')



	# -- Fourier transform to momentum space -- 
	cut=8
	CrONE = np.zeros((Ns-cut*2,4));
	NormONE = np.zeros(Ns-cut*2);
	
	## -- Two point operators ---
	for ran in range(0,Ns-cut*2):
		for i in range(cut,Ns-cut):
			for j in range(i,Ns-cut):
				if(ran == (j - i)):
					NormONE[ran] += 1
					CrONE[ran,0] += (StSt[i,j]);
					CrONE[ran,1] += (LtLt[i,j]);
					CrONE[ran,2] += (StLt[i,j]);
					CrONE[ran,3] += (JtJt[i,j]);
					#print i, j, (Sl_on[i][j]), Cr[ran][0], ran;


	for ran in range(0,Ns-cut*2):
		for i in range(0,4):
			if (NormONE[ran]==0.0):
				#print "Warning: norm = 0? at range ", ran, " and correlations ", i
				NormONE[ran] = 1.0
			CrONE[ran][i] = CrONE[ran][i]/float(NormONE[ran]);

	file1 = open("CR_pp.dat", 'w')
	for ran in range(0,Ns-cut*2):
		string = (" ".join(str(abs(CrONE[ran][i])) for i in range(0,4)))
		#print ran, string 
		outstr = str(ran)+ " " + string 
		file1.write(outstr)
		file1.write(' \n')







	### -------- plot onsite_density -------------------
	plt.figure(figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
	gs = gridspec.GridSpec(3, 1)
	gs.update(left=0.12, right=0.75, top=0.97, bottom=0.08, wspace=0.0, hspace=0.35)
	
	ax1 = plt.subplot(gs[0, :])
	plt.xlabel('site',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	titlelabel = " Lambda = "+Lambda
	plt.title(titlelabel,fontsize=16, color='black')
	plt.plot(Sx,'o-',color='black',label=r'$\langle S^x_{i} \rangle $')
	plt.plot(Sy,'o-',color='red',label=r'$\langle S^y_{i} \rangle $')
	plt.plot(Sz,'o-',color='green',label=r'$\langle S^z_{i} \rangle $')
	plt.plot((Sx+Sy+Sz)/3.0,'o-',color='magenta',label=r'$\langle M^S_{i} \rangle $')
	plt.ylim(ymax=1.05)
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	

	ax1 = plt.subplot(gs[1, :])
	plt.xlabel('site',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	titlelabel = " total-Orbital "
	#plt.title(titlelabel,fontsize=16, color='black')
	plt.plot(Lx,'o-',color='black',label=r'$\langle L^x_{i} \rangle $')
	plt.plot(Ly,'o-',color='red',label=r'$\langle L^y_{i} \rangle $')
	plt.plot(Lz,'o-',color='green',label=r'$\langle L^z_{i} \rangle $')
	plt.plot((Lx+Ly+Lz)/3.0,'o-',color='magenta',label=r'$\langle M^L_{i} \rangle $')
	#plt.ylim(ymax=1.05)
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	
	ax1 = plt.subplot(gs[2, :])
	plt.xlabel('site',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	titlelabel = " Jz "
	#plt.title(titlelabel,fontsize=16, color='black')
	plt.plot(Sz+Lz,'o-',color='black',label=r'$\langle J^z_{i} \rangle $')
	#plt.plot(Lk,'o-',color='red',label=r'$\langle S_{ib}^{z} \rangle $')
	#plt.ylim(ymax=0.6)
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	


	plt.savefig('Local.pdf')
	plt.close()
	









	### -------- plot NN and SS decay -------------------
	"""plt.figure(figsize=(8, 10), dpi=80, facecolor='w', edgecolor='k')
	gs = gridspec.GridSpec(3, 1)
	gs.update(left=0.12, right=0.75, top=0.97, bottom=0.08, wspace=0.0, hspace=0.35)"""
	plt.rc('font', family='serif')	
	fig = plt.figure(figsize=(8,6)) 
	ax=fig.add_subplot(111)

	#plt.ylabel(r'$J_{Total}/L$',size=24, labelpad=8)
	#plt.xlabel(r'$\lambda/J_{FM}$',size=24)

	markersize=6
	#ax1 = plt.subplot(gs[0, :])
	plt.xlabel('R',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	#titlelabel = " Lambda = "+Lambda
	#plt.title(titlelabel,fontsize=16, color='black')
	plt.plot(np.abs(CrONE[:,0]/CrONE[2,0]),marker="o",linewidth=2,color='black',markersize=markersize,markeredgewidth=0,label=r"$C_{S}(R)$")
	plt.plot(np.abs(CrONE[:,1]/CrONE[2,1]),marker="s",linewidth=2,color='blue',markersize=markersize,markeredgewidth=0,label=r"$C_{L}(R)$")
	#plt.plot(np.abs(CrONE[:,2]/CrONE[2,2]),'o-',color='blue',label=r'$C_{LS}(R)$')
	plt.plot(np.abs(CrONE[:,3]/CrONE[3,1]),marker="^",linewidth=2,color='red',markersize=markersize,markeredgewidth=0,label=r"$C_{J}(R)$")
	#plt.xscale('log')
	plt.yscale('log')
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


	#yticks=[0/60,10.0/60,20.0/60,30.0/60,40.0/60,50.0/60,60.0/60]
	#yticks_label=["0.0", "10/L","20/L","30/L","40/L","50/L","60/L"]
	#yticks=np.arange(-20,100,10)
	#plt.yticks(yticks,size=18)
	#plt.ylim(ymin=-0.,ymax=66)
	
	plt.xscale("log")

	plt.ylim(ymax=10)
	plt.grid(which='major',linestyle='--',alpha=0.25)
	#plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	plt.legend(loc='upper left',bbox_to_anchor=(0.74, 1.0), ncol=1, fontsize=18, frameon=False,columnspacing=1,handlelength=1.0,handletextpad=0.4)

	plt.tick_params(which=u'major',axis=u'both',length=6,color="black",direction='in')
	plt.tick_params(which=u'minor',axis=u'x',length=3,color="black",direction='in')
	plt.tick_params(which=u'minor',axis=u'y',length=3,color="black",direction='in')
	"""minor_locator = tkr.AutoMinorLocator(4)
	ax.xaxis.set_minor_locator(minor_locator)
	
	minor_locator = tkr.AutoMinorLocator(2)
	ax.yaxis.set_minor_locator(minor_locator)
	plt.grid(which='major',linestyle='--',alpha=0.25)"""

	ax.annotate("(b)", xy=(0.08,0.92), xycoords="axes fraction",fontsize=20)
	plt.savefig('NNSS_CR.pdf')
	plt.close()




	"""plt.xscale('log')
	plt.yscale('log')
	plt.ylim(ymax=10)
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	
	ax1 = plt.subplot(gs[1, :])
	plt.xlabel('site',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	titlelabel = " Real-space picture "
	#plt.title(titlelabel,fontsize=16, color='black')
	dsl = np.diagonal(StLt)
	plt.plot(dsl,'o-',color='black',label=r'$S_{i} \cdot L_{i}$')
	#plt.ylim(ymax=1,ymin=-1)
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	
	ax1 = plt.subplot(gs[2, :])
	plt.xlabel(r'$k_x$',fontsize=22)
	plt.xticks(fontsize=22)
	plt.yticks(fontsize=22)
	titlelabel = " k-space picture "
	#plt.title(titlelabel,fontsize=16, color='black')
	plt.plot(Sk/(Ns*Ns),'o-',color='black',label=r'$S(K)$')
	plt.plot(Lk/(Ns*Ns),'o-',color='red',label=r'$L(K)$')
	#plt.plot(SLk/(Ns*Ns),'o-',color='blue',label=r'$SL(K)$')
	plt.plot(Jk/(Ns*Ns),'o-',color='teal',label=r'$J(K)$')
	plt.grid(which='major',linestyle='--',alpha=0.35)
	plt.legend(loc='upper left', numpoints=1, ncol=1, fontsize=22, bbox_to_anchor=(1.00, 1.1), frameon=False,handlelength=1,handletextpad=0.4)
	plt.xticks([0.0,Ns/4,Ns/2,3*Ns/4,Ns],['0',r'$\pi/2$',r'$\pi$',r'$3\pi/2$',r'$2\pi$'])
	
	plt.savefig('NNSS_CR.pdf')
	plt.close()"""
















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
	print (string) 
	for i in range(0,num_rows):
		print (" ".join(str(x) for x in m[i]))
##### ----------------------


if __name__ == '__main__':
	sys.argv ## get the input argument
	total = len(sys.argv)
	cmdargs = sys.argv	
	main(total, cmdargs)
































