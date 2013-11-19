all :: unfilteredPDFCalc filteredPDFCalc

unfilteredPDFCalc: calcUnFilteredPDF.f90
	ifort -o unfilteredPDFCalc calcUnFilteredPDF.f90

filteredPDFCalc: calcFilteredPDF.f90
#	module load fftw/3.3.0.0
	ifort -o filteredPDFCalc calcFilteredPDF.f90 -lfftw3 -L /opt/fftw/3.3.0.0/x86_64/lib

