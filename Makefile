INC=-I /usr/lib/x86_64-linux-gnu/openmpi/include
LIB=  #-L/usr/lib/x86_64-linux-gnu -lplplotd -lplplotf95d
INC=
FC=gfortran
FCPARA=mpif90
FCOPT= -fimplicit-none -finit-local-zero -Wunused -fopenmp # -Wall -Wargument-mismatch -Wimplicit-interface
serial: Hbinitio.f90
	$(FC) $(FCOPT) Hbinitio.f90 -o Hbinitio.x -lblas -llapack
#	$(FCPARA) $(FCOPT) Hbinitio.f90 -o Hbinitio.x $(INC) -lblas -llapack -lmpi
pp: ppHbinitio.f90
	$(FC) ppHbinitio.f90 -o ppHbinitio.x $(LIB)
parallel: mpi_Hbinitio.f90
	$(FCPARA) $(FCOPT) mpi_Hbinitio.f90 -o mpi_Hbinitio.x $(INC) -lblas -llapack -lmpi
dbg: dbg.f90
	$(FCPARA) dbg.f90 -o dbg.x  $(INC) -lblas -llapack -lmpi
