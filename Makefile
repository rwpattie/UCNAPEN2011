F77    = gfortran  -Wall -Os

OBJ = penelope.o penmain.o rita.o pengeom.o timer.o penvared.o\

MOBJ = material.o \

penmain : $(OBJ)
	$(F77) -o penmain $(OBJ)
	
material : $(MOBJ)
	$(F77) -o material $(MOBJ)
	
material.o : src/material.f
	$(F77) -c src/material.f
	
penvared.o : src/penvared.f 
	$(F77) -c src/penvared.f
	
penmain.o : src/penmain.f
	$(F77) -c src/penmain.f 
	
penelope.o : src/penelope.f
	$(F77) -c src/penelope.f

rita.o : src/rita.f
	$(F77) -c src/rita.f
	
pengeom.o : src/pengeom.f
	$(F77) -c src/pengeom.f

timer.o : src/timer.f
	$(F77) -c src/timer.f

clean : 
	rm -f *.o


