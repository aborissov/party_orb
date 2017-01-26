#WARNING: this makefile compiles both RELATIVISTIC and NON-RELATIVISTIC codes
# addtional field setups are included but commented out
# CMTfields_mod contains the Guiliani et al. setup
# separatorfields_mod contains the Wilmot-Smith and Hornig (WSH) 2011 setup.
# testfields_mod contains simple test fields..
# --- make sure the right one is commented in the main body settings below

# Set the compiler flags 
FFLAGS = -O3 -J mod -ffree-line-length-0 -fbounds-check -fno-range-check -g		
#FFLAGS = -O3 -J mod -fcheck=bounds -fno-range-check -g			# original
#FFLAGS = -o0 -J mod -C -fbounds-check -g -Wall #-fdefault-real-8 	# debugging flags
#FFLAGS = -O3 -C -g -traceback -check all -warn all 			# mpif90 error flags
#FFLAGS = -O3

# name of the output executable files?
TARGETN = nrparty
TARGETR =  rparty

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

SRCDIR = src
OBJDIR = obj
BINDIR = bin
DATDIR = Data
MODDIR = mod

MODULEFLAG = -I
#FC = gfortran
#FC = gfortran44 	# oldjock
#MODULEFLAG = -module
OPFLAGS = $(QMONO) $(QSINGLE) $(QFIRSTORDER)
FC = mpif90 $(OPFLAGS)
#FC = ftn $(OPFLAGS)
PREPROFLAGS = $(NONMPIIO)


OBJFILESN = global_mod.o mpi_routines.o products_mod.o lare_functions_mod.o l3dc_fields_mod.o \
	nr_derivs_mod.o nr_rkck_mod.o nr_rkqs_mod.o nr_rkdrive_mod.o testfields_mod.o l2dc_fields_mod.o \
	l3ds_fields_mod.o sdf_common.o sdf_job_info.o sdf_control.o sdf.o sdf_output_util.o\
	sdf_input.o sdf_input_r4.o sdf_input_r8.o sdf_input_ru.o sdf_input_util.o l2ds_fields_mod.o\
	sdf_input_cartesian.o sdf_input_cartesian_r4.o sdf_input_cartesian_r8.o sdf_input_cartesian_ru.o \
	sdf_input_point.o sdf_input_point_r4.o sdf_input_point_r8.o sdf_input_point_ru.o lare_fields_mod.o\
	sdf_input_station.o sdf_input_station_r4.o sdf_input_station_r8.o sdf_input_station_ru.o\
	sdf_output.o sdf_output_r4.o sdf_output_r8.o sdf_output_ru.o sdf_output_source.o sdf_source_info.o\
	sdf_output_cartesian.o sdf_output_cartesian_r4.o sdf_output_cartesian_r8.o sdf_output_cartesian_ru.o\
	sdf_output_point.o sdf_output_point_r4.o sdf_output_point_r8.o sdf_output_point_ru.o\
	sdf_output_station.o sdf_output_station_r4.o sdf_output_station_r8.o sdf_output_station_ru.o\
	iocontrol.o input.o inputfunctions.o input_cartesian.o iocommon.o bourdinfields_mod.o\
	separatorfields_mod.o CMTfields_mod.o field_selector_mod.o gammadist_mod.o nr_main.o
OBJFILESR = global_mod.o mpi_routines.o products_mod.o lare_functions_mod.o l3dc_fields_mod.o \
	r_derivs_mod.o r_rkck_mod.o r_rkqs_mod.o r_rkdrive_mod.o r_rkcol_mod.o testfields_mod.o l2dc_fields_mod.o\
	l3ds_fields_mod.o sdf_common.o sdf_job_info.o sdf_control.o sdf.o sdf_output_util.o\
	sdf_input.o sdf_input_r4.o sdf_input_r8.o sdf_input_ru.o sdf_input_util.o l2ds_fields_mod.o\
	sdf_input_cartesian.o sdf_input_cartesian_r4.o sdf_input_cartesian_r8.o sdf_input_cartesian_ru.o \
	sdf_input_point.o sdf_input_point_r4.o sdf_input_point_r8.o sdf_input_point_ru.o lare_fields_mod.o\
	sdf_input_station.o sdf_input_station_r4.o sdf_input_station_r8.o sdf_input_station_ru.o\
	sdf_output.o sdf_output_r4.o sdf_output_r8.o sdf_output_ru.o sdf_output_source.o sdf_source_info.o\
	sdf_output_cartesian.o sdf_output_cartesian_r4.o sdf_output_cartesian_r8.o sdf_output_cartesian_ru.o\
	sdf_output_point.o sdf_output_point_r4.o sdf_output_point_r8.o sdf_output_point_ru.o\
	sdf_output_station.o sdf_output_station_r4.o sdf_output_station_r8.o sdf_output_station_ru.o\
	iocontrol.o input.o inputfunctions.o input_cartesian.o iocommon.o bourdinfields_mod.o \
	separatorfields_mod.o CMTfields_mod.o field_selector_mod.o gammadist_mod.o r_main.o


FULLTARGETN = $(BINDIR)/$(TARGETN)
FULLTARGETR = $(BINDIR)/$(TARGETR)

VPATH = $(SRCDIR):$(OBJDIR):$(SRCDIR)/core:$(SRCDIR)/core/othermods:$(SRCDIR)/core/nr_rkmods:\
	$(SRCDIR)/core/r_rkmods:$(SRCDIR)/core/laremods:$(SRCDIR)/core/laremods/cfd:$(SRCDIR)/core/laremods/sdf

DATDIRN = $(DATDIR)/DataN
DATDIRR = $(DATDIR)/DataR

# Rule to build the fortran files
%.o: %.f90
	@mkdir -p $(BINDIR) $(OBJDIR) $(DATDIR) $(DATDIRN) $(DATDIRR) $(MODDIR)
	$(FC) -c $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<	
%.o: %.F90
	@mkdir -p $(BINDIR) $(OBJDIR) 
	$(FC) -c $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $(OBJDIR)/$@ $<

all: $(FULLTARGETN) $(FULLTARGETR)

$(FULLTARGETN): $(OBJFILESN)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILESN))

$(FULLTARGETR): $(OBJFILESR)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILESR))

nonrel: .$(OBJFILESN)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILESN))

rel:  $(OBJFILESR)
	$(FC) $(FFLAGS) $(MODULEFLAG) $(OBJDIR) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILESR))

clean:
	@rm -rf *~ $(BINDIR) $(OBJDIR) $(MODDIR) $(SRCDIR)/*~ $(SRCDIR)/core/*~

tidy:
	@rm -rf $(OBJDIR) $(SRCDIR)/*~

datatidy:
	@rm -rf $(DATDIR)/*

.PHONY: clean tidy datatidy rel nonrel

# All the dependencies
global_mod.o: global_mod.f90 sdf_job_info.o
products_mod.o: products_mod.f90 global_mod.o
gammadist_mod.o: gammadist_mod.f90 global_mod.o
lare_functions_mod.o: lare_functions_mod.f90 global_mod.o
iocontrol.o: iocontrol.f90 global_mod.o iocommon.o input.o 
input.o: input.f90 global_mod.o iocommon.o inputfunctions.o 
inputfunctions.o: inputfunctions.f90 global_mod.o iocommon.o  
iocommon.o: iocommon.f90 global_mod.o 
input_cartesian.o: input_cartesian.f90 iocommon.o inputfunctions.o 
testfields_mod.o: testfields_mod.f90 global_mod.o products_mod.o
separatorfields_mod.o: separatorfields_mod.f90 global_mod.o products_mod.o
CMTfields_mod.o: CMTfields_mod.f90 global_mod.o products_mod.o
l3dc_fields_mod.o: l3dc_fields_mod.f90 lare_functions_mod.o global_mod.o iocommon.o iocontrol.o input.o input_cartesian.o #separatorfields_mod.o
l2dc_fields_mod.o: l2dc_fields_mod.f90 lare_functions_mod.o global_mod.o iocommon.o iocontrol.o input.o input_cartesian.o #separatorfields_mod.o
lare_fields_mod.o: lare_fields_mod.f90 lare_functions_mod.o global_mod.o
sdf_common.o: sdf_common.f90 sdf_job_info.o global_mod.o
sdf_job_info.o: sdf_job_info.f90
sdf_source_info.o: sdf_source_info_dummy.f90
sdf_input.o: sdf_input.f90 sdf_input_r4.o sdf_input_r8.o
sdf_input_r4.o: sdf_input_r4.f90 sdf_input_ru.o
sdf_input_r8.o: sdf_input_r8.f90 sdf_input_ru.o
sdf_input_ru.o: sdf_input_ru.f90 sdf_common.o
sdf_input_cartesian.o: sdf_input.f90 sdf_input_cartesian_r4.o sdf_input_cartesian_r8.o
sdf_input_cartesian_r4.o: sdf_input_cartesian_r4.f90 sdf_input_cartesian_ru.o
sdf_input_cartesian_r8.o: sdf_input_cartesian_r8.f90 sdf_input_cartesian_ru.o
sdf_input_cartesian_ru.o: sdf_input_cartesian_ru.f90 sdf_input_ru.o
sdf_input_point.o: sdf_input.f90 sdf_input_point_r4.o sdf_input_point_r8.o
sdf_input_point_r4.o: sdf_input_point_r4.f90 sdf_input_point_ru.o
sdf_input_point_r8.o: sdf_input_point_r8.f90 sdf_input_point_ru.o
sdf_input_point_ru.o: sdf_input_point_ru.f90 sdf_input_ru.o
sdf_input_station.o: sdf_input.f90 sdf_input_station_r4.o sdf_input_station_r8.o
sdf_input_station_r4.o: sdf_input_station_r4.f90 sdf_input_station_ru.o
sdf_input_station_r8.o: sdf_input_station_r8.f90 sdf_input_station_ru.o
sdf_input_station_ru.o: sdf_input_station_ru.f90 sdf_input_ru.o
sdf_input_util.o: sdf_input_util.f90 sdf_input.o sdf_input_cartesian.o sdf_input_point.o sdf_input_station.o sdf_output_station_ru.o
sdf_output.o: sdf_output.f90 sdf_output_r4.o sdf_output_r8.o sdf_output_ru.o
sdf_output_r4.o: sdf_output_r4.f90 sdf_output_ru.o
sdf_output_r8.o: sdf_output_r8.f90 sdf_output_ru.o
sdf_output_ru.o: sdf_output_ru.f90 sdf_common.o
sdf_output_cartesian.o: sdf_output.f90 sdf_output_cartesian_r4.o sdf_output_cartesian_r8.o
sdf_output_cartesian_r4.o: sdf_output_cartesian_r4.f90 sdf_output_cartesian_ru.o
sdf_output_cartesian_r8.o: sdf_output_cartesian_r8.f90 sdf_output_cartesian_ru.o
sdf_output_cartesian_ru.o: sdf_output_cartesian_ru.f90 sdf_output_ru.o
sdf_output_point.o: sdf_output.f90 sdf_output_point_r4.o sdf_output_point_r8.o
sdf_output_point_r4.o: sdf_output_point_r4.f90 sdf_output_point_ru.o
sdf_output_point_r8.o: sdf_output_point_r8.f90 sdf_output_point_ru.o
sdf_output_point_ru.o: sdf_output_point_ru.f90 sdf_output_ru.o sdf_common.o
sdf_output_station.o: sdf_output_station.f90 sdf_output_station_r4.o sdf_output_station_r8.o
sdf_output_station_r4.o: sdf_output_station_r4.f90 sdf_output_station_ru.o
sdf_output_station_r8.o: sdf_output_station_r8.f90 sdf_output_station_ru.o
sdf_output_station_ru.o: sdf_output_station_ru.f90 sdf_output_ru.o
sdf_output_util.o: sdf_output_util.f90 sdf_output_cartesian_ru.o sdf_output_point_ru.o 
sdf_output_source.o: sdf_output_source.f90 sdf_common.o sdf_source_info.o sdf_output.o
sdf_control.o: sdf_control.f90 sdf_output_util.o
sdf.o: sdf.f90 sdf_control.o sdf_input.o sdf_input_cartesian.o sdf_input_point.o sdf_input_station.o sdf_input_util.o\
sdf_output.o sdf_output_cartesian.o sdf_output_point.o sdf_output_station.o sdf_output_source.o sdf_md5.o
mpi_routines.o:mpi_routines.f90 global_mod.o l3dc_fields_mod.o l3ds_fields_mod.o l2dc_fields_mod.o l2ds_fields_mod.o
l3ds_fields_mod.o: l3ds_fields_mod.f90 sdf_common.o global_mod.o sdf.o sdf_job_info.o lare_functions_mod.o
l2ds_fields_mod.o: l2ds_fields_mod.f90 sdf_common.o global_mod.o sdf.o sdf_job_info.o lare_functions_mod.o
bourdinfields_mod.o: bourdinfields_mod.f90 lare_functions_mod.o global_mod.o products_mod.o
field_selector_mod.o: field_selector_mod.f90 lare_fields_mod.o lare_functions_mod.o CMTfields_mod.o separatorfields_mod.o testfields_mod.o bourdinfields_mod.o global_mod.o
##--NREL DEPENDENCIES
nr_derivs_mod.o: nr_derivs_mod.f90 global_mod.o field_selector_mod.o products_mod.o 
nr_rkck_mod.o: nr_rkck_mod.f90 nr_derivs_mod.o global_mod.o field_selector_mod.o
nr_rkqs_mod.o: nr_rkqs_mod.f90 global_mod.o nr_rkck_mod.o field_selector_mod.o
nr_rkdrive_mod.o: nr_rkdrive_mod.f90 global_mod.o nr_derivs_mod.o nr_rkqs_mod.o field_selector_mod.o
#--REL DEPENDENCIES
r_rkcol_mod.o: r_rkcol_mod.f90 global_mod.o products_mod.o r_derivs_mod.o 
r_derivs_mod.o: r_derivs_mod.f90 global_mod.o products_mod.o field_selector_mod.o
r_rkck_mod.o: r_rkck_mod.f90 r_derivs_mod.o global_mod.o field_selector_mod.o
r_rkqs_mod.o: r_rkqs_mod.f90 global_mod.o r_rkck_mod.o r_rkcol_mod.o field_selector_mod.o
r_rkdrive_mod.o: r_rkdrive_mod.f90 global_mod.o r_derivs_mod.o r_rkqs_mod.o field_selector_mod.o
#mp
nr_main.o: nr_main.f90 global_mod.o mpi_routines.o nr_rkdrive_mod.o products_mod.o field_selector_mod.o lare_functions_mod.o bourdinfields_mod.o gammadist_mod.o
r_main.o: r_main.f90 global_mod.o mpi_routines.o r_rkdrive_mod.o products_mod.o field_selector_mod.o lare_functions_mod.o bourdinfields_mod.o gammadist_mod.o
