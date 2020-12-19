#
# Makefile for caffa3d program
#

# Compiler flags:
LFLAGS = -llapack
#debug
#F90FLAGS = -Wall -g -fstack-protector-all
#realse
F90FLAGS = -Wall -g
# Compiler:
F90 = gfortran

MOD_FILES=\
    m_allocatable.f90\
    m_utils.f90\
    m_matrix.f90\
    m_geometry.f90\
    m_output.f90\
    m_fieldInitialization.f90\
    m_sparseMatrix.f90\
    m_gradients.f90\
    m_interpolation.f90\
    m_fieldManipulation.f90\
    m_scalarFluxes.f90\
    m_temperature.f90\
    m_concentration.f90\
    m_mhd.f90\
    m_facefluxVelocity.f90\

LINEAR_SOLVER_FILES=\
    iccg.f90 \
    bicgstab.f90 \
    dpcg.f90 \
    mgmres.f90 \
    gauss-seidel.f90\
	linear_solvers.f90

TURBULENCE_FILES=\
    k_epsilon_std.f90\
    k_epsilon_std_2lewt.f90\
    k_omega_SST.f90\
    k_omega_SST_DDES.f90\
    k_omega_SST_IDDES.f90\
    k_epsilon_rng.f90\
    k_epsilon_rlzb.f90\
    k_epsilon_rlzb_2lewt.f90\
    k_eqn_eddy.f90\
    m_spalart_allmaras.f90\
    m_smagorinsky.f90

SRC_FILES=\
    asm_heatflux_terms.f90\
    adjustMassFlow.f90\
    bcin.f90\
    bpres.f90\
    faceflux_mass.f90\
    calcheatflux.f90\
    continuityErrors.f90\
    calcp_simple.f90\
    calcp_piso.f90\
    calcstress.f90\
    calc_strain_and_vorticity.f90\
    calcuvw.f90\
    calc_statistics.f90\
    CourantNo.f90\
    constant_mass_flow_forcing.f90\
    create_fields.f90\
    updateVelocityAtBoundary.f90\
    modify_viscosity.f90\
    modify_viscosity_inlet.f90\
    laplacian.f90\
    init.f90\
    readfiles.f90\
    read_input.f90\
    time_shift.f90\
    writefiles.f90\
    write_restart_files.f90\
    writehistory.f90\
    main.f90 

#
# How to create object files:
# 
MODS = ${MOD_FILES:.f90=.o}
TURBULENCE = ${TURBULENCE_FILES:.f90=.o}
LINEAR_SOLVERS = ${LINEAR_SOLVER_FILES:.f90=.o}
SRCOBJS = ${SRC_FILES:.f90=.o}

##################################################################
# Targets for make.
##################################################################

all: cappuccino

cappuccino: ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${SRCOBJS}
	@echo  "Linking" $@ "... "
	${F90} ${SRCOBJS} ${MODS} ${TURBULENCE} ${LINEAR_SOLVERS} ${LFLAGS} ${INCS} -o ../bin/cappuccino 

.PHONY: clean
clean:
	@rm -f *.o *.mod ../bin/cappuccino

##################################################################
# Generic rules
##################################################################

.SUFFIXES : .f90

.f90.o:
	${F90} ${F90FLAGS} -c ${INCS}  ${@:.o=.f90}
