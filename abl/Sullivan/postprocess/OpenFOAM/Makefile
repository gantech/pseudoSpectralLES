all :openFoamBCs openFoamSideInletBCs openFoamSideOutletBCs openFoamTopBCs openFoamICs 

openFoamBCs : writeOpenFOAMBCs.f90
	ifort -o openFoamBCs writeOpenFOAMBCs.f90

openFoamSideInletBCs : writeOpenFOAMsideInletBCs.f90
	ifort -o openFoamSideInletBCs writeOpenFOAMsideInletBCs.f90

openFoamSideOutletBCs : writeOpenFOAMsideOutletBCs.f90
	ifort -o openFoamSideOutletBCs writeOpenFOAMsideOutletBCs.f90

openFoamTopBCs : writeOpenFOAMtopBCs.f90
	ifort -o openFoamTopBCs writeOpenFOAMtopBCs.f90

openFoamICs : writeOpenFOAMInitialConditions.f90
	ifort -o openFoamICs writeOpenFOAMInitialConditions.f90
