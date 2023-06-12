F08 = gfortran
FFLAGS = -O3 -g -Wall -Warray-bounds -ffixed-line-length-none -fbounds-check -fopenmp -O3 -ffast-math -funroll-loops -march=native -mtune=native -J mod
SRCDIR = ./src/pntwrks
SLVRDIR = ./src/slvrs
MAINDIR = ./src
MODDIR = ./mod
OBJDIR = ./obj

#####################################################################
OBJS =\
         structs.o pntst.o eqslvrs.o krnl_cmn.o krnl_sph.o krnl_wls.o krnl_rbf.o krnl_mls.o krnl_krg.o krnl_gfd.o bc.o intrf.o trnsprt.o ns.o ns_cbs.o sdf.o prmtrs.o slv_trnsprt.o slv_trnsprt_ss.o slv_ch.o slv_ac.o slv_ls.o slv_ns_1phs.o slv_ns_2phs.o slv_ns_lgr.o slv_elst.o slv_intrp.o slv_sdf.o slv_stfn.o main.o

#####################################################################

RUN: $(OBJS)
	$(F08) $(FFLAGS) $(OBJS) -o run

clean:
	rm *.mod *.o *~ run

#####################################################################
structs.o        	       : $(SRCDIR)/structs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/structs.f08

pntst.o        	       : $(SRCDIR)/pntst.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/pntst.f08

eqslvrs.o        	       : $(SRCDIR)/eqslvrs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/eqslvrs.f08

krnl_cmn.o        	       : $(SRCDIR)/krnl_cmn.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_cmn.f08

krnl_sph.o        	       : $(SRCDIR)/krnl_sph.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_sph.f08

krnl_wls.o        	       : $(SRCDIR)/krnl_wls.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_wls.f08

krnl_rbf.o        	       : $(SRCDIR)/krnl_rbf.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_rbf.f08

krnl_mls.o        	       : $(SRCDIR)/krnl_mls.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_mls.f08

krnl_krg.o        	       : $(SRCDIR)/krnl_krg.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_krg.f08

krnl_gfd.o        	       : $(SRCDIR)/krnl_gfd.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/krnl_gfd.f08

bc.o        	       : $(SRCDIR)/bc.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/bc.f08

intrf.o        	       : $(SRCDIR)/intrf.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/intrf.f08

trnsprt.o        	       : $(SRCDIR)/trnsprt.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/trnsprt.f08

ns.o        	       : $(SRCDIR)/ns.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/ns.f08

ns_cbs.o        	       : $(SRCDIR)/ns_cbs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/ns_cbs.f08

sdf.o        	       : $(SRCDIR)/sdf.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/sdf.f08

prmtrs.o        	       : $(SRCDIR)/prmtrs.f08
	$(F08) $(FFLAGS)  -c $(SRCDIR)/prmtrs.f08
#####################################################################
slv_trnsprt.o        	       : $(SLVRDIR)/slv_trnsprt.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_trnsprt.f08

slv_trnsprt_ss.o        	       : $(SLVRDIR)/slv_trnsprt_ss.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_trnsprt_ss.f08

slv_ls.o        	       : $(SLVRDIR)/slv_ls.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ls.f08

slv_ch.o        	       : $(SLVRDIR)/slv_ch.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ch.f08

slv_ac.o        	       : $(SLVRDIR)/slv_ac.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ac.f08

slv_ns_1phs.o        	       : $(SLVRDIR)/slv_ns_1phs.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ns_1phs.f08

slv_ns_2phs.o        	       : $(SLVRDIR)/slv_ns_2phs.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ns_2phs.f08

slv_ns_lgr.o        	       : $(SLVRDIR)/slv_ns_lgr.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_ns_lgr.f08

slv_elst.o        	       : $(SLVRDIR)/slv_elst.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_elst.f08

slv_intrp.o        	       : $(SLVRDIR)/slv_intrp.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_intrp.f08

slv_sdf.o        	       : $(SLVRDIR)/slv_sdf.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_sdf.f08
	
slv_stfn.o        	       : $(SLVRDIR)/slv_stfn.f08
	$(F08) $(FFLAGS)  -c $(SLVRDIR)/slv_stfn.f08
#####################################################################

main.o        	       : $(MAINDIR)/main.f08
	$(F08) $(FFLAGS)  -c $(MAINDIR)/main.f08

#####################################################################
