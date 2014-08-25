all :openFoamBCs openFoamSideInletBCs openFoamSideOutletBCs openFoamICs 

openFoamBCs : writeOpenFOAMBCs.f90
	ifort -o openFoamBCs writeOpenFOAMBCs.f90

openFoamSideInletBCs : writeOpenFOAMsideInletBCs.f90
	ifort -o openFoamSideInletBCs writeOpenFOAMsideInletBCs.f90

openFoamSideOutletBCs : writeOpenFOAMsideOutletBCs.f90
	ifort -o openFoamSideOutletBCs writeOpenFOAMsideOutletBCs.f90

openFoamICs : writeOpenFOAMInitialConditions.f90
	ifort -o openFoamICs writeOpenFOAMInitialConditions.f90
