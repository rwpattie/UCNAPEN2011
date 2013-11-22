F77    = gfortran -Wall -O3 -fno-align-commons

CLIBS ='/usr/lib/i386-linux-gnu'

ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)   -lNew -lMinuit -lSpectrum
ROOTGLIBS     = $(shell root-config --glibs)
CXX           = gcc
CXXFLAGS      = -O3 -Wall -fPIC -g -Wl,--no-as-needed
SOFLAGS       = -shared
CXXFLAGS     += $(ROOTCFLAGS)
LIBS          = $(ROOTLIBS)
GLIBS         = $(ROOTGLIBS)
GLIBS         +=  -lz -lRMySQL -lnsl -lcrypt -ldl -lstdc++

%.o: src/%.f
	$(F77) -c $< $(GLIBS)
	
root_util_io.o: src/root_util_io.C
	$(CXX) $(CXXFLAGS) -c $< $(GLIBS)
	 
OBJ = penelope.o ucnapenmain.o rita.o pengeom.o timer.o penvared.o penfield.o b_field.o\
      pen_math_utils.o pen_io_utils.o ucna_detector_utils.o  ucna_paw_io.o \
      full_neutron.o al_decay.o tin_decay.o bismuth.o cd_decay.o\
      gamma_source.o cu_capture.o in_decay.o root_util_io.o\

MOBJ = material.o \

GOBJ = geo_writer.o \
	
all : geo_writer ucnapenmain material

ucnapenmain : $(OBJ)
	$(F77) -o ucnapenmain $(OBJ) \
	$(CLIBS)/libpacklib.so \
	$(CLIBS)/libpawlib.so $(GLIBS)
	
geo_writer : src/geo_writer.f
	$(F77) src/geo_writer.f -o geo_writer

material : $(MOBJ)
	$(F77) -o material $(MOBJ)

clean : 
	rm -f *.o
	rm ucnapenmain
	rm material
	rm geo_writer


