from pylab import loadtxt, transpose, size, logspace, ones
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
    import matplotlib.pyplot as p
    
    with open('wSpectrumPeakWavenumber','w') as wPeakFile:
        wPeakFile.write('iz           wSpectrumPeakWavenumber'+'\n')
        for z in range(size(zLevels)):
            wPeakFile.write(str(z+1)+'           '+str(w[z,:].argmax())+'\n')


    print logspace(-7,3,11)
    #Make spectra plots at each level
    for z in range(size(zLevels)):
        p.figure()
        p.loglog(wavenumber,u[z,:],label='U')
        p.loglog(wavenumber,v[z,:],label='V')
        p.loglog(wavenumber,w[z,:],label='W')
        km5by3 = 200.0*wavenumber**(-5.0/3.0) 
        p.loglog(wavenumber, km5by3,'k--',label='km5by3')
        p.loglog(w[z,:].argmax() * ones(11), logspace(-7,3,11),'k')
#        p.xlim(0,0.01)
        p.title('z = '+str((z+0.5)*dz)+'m')
        p.xlabel('k (in 1/m')
        p.ylabel('Amplitude')
        p.savefig('z'+str(z)+'_Spectra.png')
    
        
