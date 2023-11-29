#!/bin/sh -f

FC=ifort
FFLAGS=-132 -r8 -O3 -check bounds
LIB_NetCDF=-L$(NETCDF_LIB) -lnetcdff -lnetcdf
INC_NetCDF=-I$(NETCDF_INC)

FMODS=
FOBJS=src2nc.o 

.SUFFIXES:.o .f .f90 .F90
.f.o:
	$(FC) -c $(FFLAGS) $<
.f90.o:
	$(FC) -c $(FFLAGS) $(INC_NetCDF) $<
#
#----------------------------------------------
#  Tasks
#----------------------------------------------
#
all:src2nc.x 
	@echo '---------------------------->>> END OF COMPILATION'

src2nc.x:src2nc.o
	$(FC) $(FFLAGS) $^ -o $@ $(LIB_NetCDF)

new:
	@rm -rf *.o *.mod core*
	@make

clean:
	@rm -f core core.* *~ *.o *.mod

