from pylab import loadtxt, transpose, size
nLevels = 50
nx = 768
xl = 5120.0
data = loadtxt('u_spectra.dat')
wavenumber = data[:,0]#/xl
u = transpose(data[:,1:])
data = loadtxt('v_spectra.dat')
v = transpose(data[:,1:])
data = loadtxt('w_spectra.dat')
w = transpose(data[:,1:])
zLevels = range(1,51,1)
dz = 8.0
del data

if __name__=="__main__":
    import matplotlib
    import matplotlib.pyplot as p
    matplotlib.use('Agg')
    
    #Make spectra plots at each level
    for z in range(size(zLevels)):
        p.figure()
        p.semilogx(wavenumber,u[z,:],label='U')
        p.semilogx(wavenumber,v[z,:],label='V')
        p.semilogx(wavenumber,w[z,:],label='W')
#        p.xlim(0,0.01)
        p.title('z = '+str((z+0.5)*dz)+'m')
        p.xlabel('k (in 1/m')
        p.ylabel('Amplitude')
        p.savefig('z'+str(z)+'_Spectra.png')
        
