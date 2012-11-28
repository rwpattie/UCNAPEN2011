F77    = gfortran  -Wall -O3
CLIBS ='/usr/lib'
OBJ = penelope.o ucnapenmain.o rita.o pengeom.o timer.o penvared.o penfield.o b_field.o\
      pen_math_utils.o pen_io_utils.o ucna_detector_utils.o ucna_paw_io.o\
      full_neutron.o al_decay.o tin_decay.o bismuth.o cd_decay.o \
      gamma_source.o cu_capture.o\

MOBJ = material.o \

tester : tester.o
	$(F77) -o tester tester.o \
	$(CLIBS)/libpacklib.so \
	$(CLIBS)/libpawlib.so -lnsl -lcrypt -ldl

ucnapenmain : $(OBJ)
	$(F77) -o ucnapenmain $(OBJ) \
	$(CLIBS)/libpacklib.so \
	$(CLIBS)/libpawlib.so -lnsl -lcrypt -ldl
	
material : $(MOBJ)
	$(F77) -o material $(MOBJ)

tester.o : src/tester.f
	$(F77) -c src/tester.f
	
material.o : src/material.f
	$(F77) -c src/material.f
	
ucna_detector_utils.o : src/ucna_detector_utils.f
	$(F77) -c src/ucna_detector_utils.f
	
ucna_paw_io.o : src/ucna_paw_io.f
	$(F77) -c src/ucna_paw_io.f
	
pen_math_utils.o : src/pen_math_utils.f
	$(F77) -c src/pen_math_utils.f
	
pen_io_utils.o : src/pen_io_utils.f
	$(F77) -c src/pen_io_utils.f
	
penfield.o : src/penfield.f
	$(F77) -c src/penfield.f
	
penvared.o : src/penvared.f 
	$(F77) -c src/penvared.f
	
ucnapenmain.o : src/ucnapenmain.f
	$(F77) -c src/ucnapenmain.f 
	
penelope.o : src/penelope.f
	$(F77) -c src/penelope.f

rita.o : src/rita.f
	$(F77) -c src/rita.f

b_field.o : src/b_field.f
	$(F77) -c src/b_field.f
	
pengeom.o : src/pengeom.f
	$(F77) -c src/pengeom.f

timer.o : src/timer.f
	$(F77) -c src/timer.f
	
al_decay.o : src/al_decay.f
	$(F77) -c src/al_decay.f
	
cd_decay.o : src/cd_decay.f
	$(F77) -c src/cd_decay.f
	
bismuth.o : src/bismuth.f
	$(F77) -c src/bismuth.f

tin_decay.o : src/tin_decay.f
	$(F77) -c src/tin_decay.f

full_neutron.o : src/full_neutron.f
	$(F77) -c src/full_neutron.f
	
cu_capture.o : src/cu_capture.f
	$(F77) -c src/cu_capture.f

gamma_source.o : src/gamma_source.f
	$(F77) -c src/gamma_source.f



clean : 
	rm -f *.o


