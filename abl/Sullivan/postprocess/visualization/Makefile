all: isocontoursFiltered isocontoursUnfiltered

isocontoursFiltered : generateIsocontoursFiltered.f90 LIB_VTK_IO.f90
#	module load fftw3
	ifort -I${TACC_FFTW3_INC} -o isocontoursFiltered generateIsocontoursFiltered.f90 LIB_VTK_IO.f90 -lfftw3 -L${TACC_FFTW3_LIB}

isocontoursUnfiltered : generateIsocontoursUnfiltered.f90 LIB_VTK_IO.f90
#	module load fftw3
	ifort -I${TACC_FFTW3_INC} -o isocontoursUnfiltered generateIsocontoursUnfiltered.f90 LIB_VTK_IO.f90 -L${TACC_FFTW3_LIB}

hdf5files: generateHDF5Unfiltered.f90 
	ifort -I${TACC_HDF5_INC} -o hdf5Unfiltered generateHDF5Unfiltered.f90 -L${TACC_HDF5_LIB} -lhdf5_fortran
