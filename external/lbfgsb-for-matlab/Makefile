# Linux settings.
MATLAB_HOME = /usr/lib/matlab-8.3/
MEX         = $(MATLAB_HOME)/bin/mex
MEXSUFFIX   = mexa64
CXX         = g++-4.7
F77         = gfortran-4.7
CFLAGS      = -O3 -fPIC -pthread 
FFLAGS      = -O3 -fPIC -fexceptions 

TARGET = lbfgsb.$(MEXSUFFIX)
OBJS   = solver.o matlabexception.o matlabscalar.o matlabstring.o   \
         matlabmatrix.o arrayofmatrices.o program.o matlabprogram.o \
         lbfgsb.o

CFLAGS += -Wall -ansi -DMATLAB_MEXFILE

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(MATLAB_HOME)/extern/include -o $@ -c $^

%.o: %.f
	$(F77) $(FFLAGS) -o $@ -c $^

$(TARGET): $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) FC=$(FCC) LD=$(CXX) -lgfortran -lm \
        -O -output $@ $^

clean:
	rm -f *.o $(TARGET)

