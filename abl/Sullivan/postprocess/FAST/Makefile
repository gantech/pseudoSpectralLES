all: FASTfilesWrite FASTfilesWriteCWF1p0 

FASTfilesWrite: writeFASTfiles.f90
	ftn -o FASTfilesWrite writeFASTfiles.f90
#	ifort -o FASTfilesWrite writeFASTfiles.f90 # Use sizeOfReal=1 in intel. Reads in words.
#	pgf90 -o FASTfilesWrite writeFASTfiles.f90 # Use sizeOfReal=4 in pgi. Reads in bytes.

FASTfilesWriteCWF1p0: writeFASTfiles_CWF1p0.f90
	ftn -o FASTfilesWriteCWF1p0 writeFASTfiles_CWF1p0.f90

