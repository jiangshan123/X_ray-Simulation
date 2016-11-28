# variables to note the compiler and options
FCOMP := g++
NODEBUG := 

# set the compiler flags depending on whether we are in DEBUG mode
ifdef NODEBUG
  FFLAGS := -O
else
  FFLAGS := -g
endif

# the list of source files that need to be compiled to make the executable
c++sources += 3D_math.cpp X-Ray_Simulator.cpp main.cpp 

# the list of object files that need to be created
c++objects = $(c++sources:.cpp=.o)

# default target
all: X-Ray_Simulator

# inter-file dependencies
include GMake.inc

# rule to make a .o from a .cpp and .h file
%.o:%.cpp %.h
	$(FCOMP) $(FFLAGS) -c $<

# rule to build the main target
X-Ray_Simulator: $(c++objects)
	$(FCOMP) $(FFLAGS) $(c++objects) -o X-Ray_Simulator 

# a special target that cleans up
clean:
	rm -f *.o *.mod

# PHONY tells make that clean is not a real file, so it will ensure
# that if you type "make clean" this is always executed.  This is only
# a problem if, somehow, you create an actualy file named "clean" in
# your directory -- then make will think that clean is up to date.
.PHONY: clean
