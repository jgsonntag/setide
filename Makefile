#
#  Makefile for setide program
#
CFLAGS = -lm
SRC = setide.cpp 
OBJ = setide.o
ULIBS = /home/sonntag/Libcpp/libjohn2.a


/home/sonntag/bin/setide : setide
	mv setide /home/sonntag/bin/setide

setide: $(OBJ) $(ULIBS)
	g++ $(CFLAGS) -L/home/sonntag/Libcpp -o setide $(OBJ) -ljohn2
	
setide.o: setide.cpp
	g++ -c setide.cpp

$(ULIBS): FORCE
	cd /home/sonntag/Libcpp; $(MAKE)

FORCE:
