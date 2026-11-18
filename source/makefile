# Makefile for ecomsed 

.SUFFIXES: .a .o .f

FC = ifort
LIBS = -L/usr/local/lib -lnetcdf -lnetcdff -lhdf5 -lz -L/usr/lib64
FFLAGS = -mcmodel=medium -fp-model precise -O3 -pg -I/usr/local/include
LFLAGS = -mcmodel=medium -fp-model precise -O3 -pg -I/usr/local/include

OBJECTS=  bcdata.o first.o disply.o maxmin.o print.o prints.o prtxy.o \
          setdom.o sinter.o slicexz.o sliceyz.o tands.o zeroes.o reordr.o \
          sedflx.o sedic.o suslod.o fhcalc.o wavesmb.o wavedon.o \
          chembed.o randb.o oclamp.o \
          advave.o advq.o advt.o \
          antidif.o antidif_chem.o antidif_sed.o \
          archive.o baropg.o pclamp.o ptide_noaa.o\
          alphabc.o densityUNESCO.o ecom3d.o \
          profu.o profv.o proft.o profq.o  profe.o\
          smag.o transport.o vertvl.o advsed.o profsed.o bulk.o fzol.o jday.o \
          longwave.o ncld.o n_cloud.o qsat.o vapor.o wreal.o traninp.o \
          chemflx.o partrak.o onepart.o advcon.o stress.o bottau.o advchem.o \
          datehr.o getdate.o hitflx.o \
          north.o solar.o zenith.o mixed.o

OBJECTS2= putcdf.o tscdf.o \
	extrnl.o advu.o advv.o bcond.o setdom_corners.o epcpercell.o \
   smooth.o wandd_kky.o read_corner_loc.o update_cbc.o wetbulb.o airdensity.o update_cbc_mannings.o \
	 netcdf_read_restart.o netcdf_read_synop_met.o netcdf_read_ifric2dt.o \
	 netcdf_read_corner_loc.o netcdf_read_model_grid.o netcdf_read_bc_elev.o \
	 netcdf_read_bc_ts.o netcdf_read_bc_diff1.o netcdf_read_bc_diff2.o \
	 netcdf_read_ic_ts.o netcdf_read_bc_met.o netcdf_read_bc_dis.o


s--pom.x: ${OBJECTS} $(OBJECTS2)
	$(FC) $(LFLAGS) -o s--pom.x $(OBJECTS) $(OBJECTS2) $(LIBS)
	@echo "new s--pom.x model is now up to date."

$(OBJECTS): comdeck makefile
$(OBJECTS2): comdeck netcdf.inc makefile

# Cleaning target
clean:
	@rm -f *.o *.oo
