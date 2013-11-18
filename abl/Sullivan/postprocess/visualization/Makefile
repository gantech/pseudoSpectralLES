all: isocontoursFiltered isocontoursUnfiltered

isocontoursFiltered : generateIsocontoursFiltered.f90 LIB_VTK_IO.f90
	ifort -o isocontoursFiltered generateIsocontoursFiltered.f90 LIB_VTK_IO.f90 -lfftw3 -L /opt/fftw/3.3.0.0/x86_64/lib

isocontoursUnfiltered : generateIsocontoursUnfiltered.f90 LIB_VTK_IO.f90
	ifort -o isocontoursUnfiltered generateIsocontoursUnfiltered.f90 LIB_VTK_IO.f90
