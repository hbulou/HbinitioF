INC=-I /usr/lib/x86_64-linux-gnu/openmpi/include
INC=
FC=gfortran
FCPARA=mpif90
FCOPT= -fimplicit-none # -Wall -Wargument-mismatch -Wimplicit-interface
serial: Hbinitio.f90
	$(FC) Hbinitio.f90 -o Hbinitio.x -lblas -llapack
parallel: mpi_Hbinitio.f90
	$(FCPARA) $(FCOPT) mpi_Hbinitio.f90 -o mpi_Hbinitio.x $(INC) -lblas -llapack -lmpi
dbg: dbg.f90
	$(FCPARA) dbg.f90 -o dbg.x  $(INC) -lblas -llapack -lmpi
