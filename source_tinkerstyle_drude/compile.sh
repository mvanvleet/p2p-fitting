#!/bin/bash

#OPT="-check bounds -check uninit -check format -traceback "
OPT=""
#nr="nrtype.f90 nrutil.f90 nr.f90 mnbrak.f90 brent.f90 linmin.f90 frprmn.f90"
nr="nrtype.f90 nrutil.f90 nr.f90 svd_routines.f90"

ifort $OPT $nr variable.f90 electrostatic_interaction_functions.f90 routines.f90 scf_drude.f90 amoeba.f90 functions.f90  main_p2p.f90 -o main