# Linux settings.
MATLAB_HOME = /usr/lib/matlab-8.4/
MEX         = $(MATLAB_HOME)/bin/mex
MEXSUFFIX   = mexa64
CXX         = g++-4.7
CFLAGS      = -O3 -fPIC -pthread 

TARGET = ErrorEvaluation.$(MEXSUFFIX)
OBJS   = ErrorEvaluation.o GeneralizedProcrustes.o \
	patternRecognitionPCA.o Mle.o \
	GaussVector.o UnsupervisedLearning.o

CFLAGS += -Wall -ansi -DMATLAB_MEXFILE

all: $(TARGET)

%.o: %.cpp
	$(CXX) $(CFLAGS) -I$(MATLAB_HOME)/extern/include  -o $@ -c $^

$(TARGET): $(OBJS)
	$(MEX) -cxx CXX=$(CXX) CC=$(CXX) LD=$(CXX) $(MATLAB_HOME)/bin/glnxa64/libmwlapack.so $(MATLAB_HOME)/bin/glnxa64/libmwblas.so \
        -O -output $@ $^

clean:
	rm -f *.o $(TARGET)
