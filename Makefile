CXX = g++
CC  = gcc
FF  = gfortran
Opt = -Ofast

include Makefile.in

CXX_FLAGS = -ansi -fexceptions -DMATLAB_MEX_FILE -std=c++11 -fopenmp -march=native \
			-D_GNU_SOURCE -fPIC -fno-omit-frame-pointer -Wno-write-strings -pthread\
			$(Opt) -DNDEBUG -fopenmp -ffast-math
			
CXX_INCLUDE = -I./include/mexplus \
              -I./include/ \
              -I$(MATLAB_ROOT)extern/include \
              -I$(MATLAB_ROOT)simulink/include
              
MATLAB_LINKS = $(Opt) -pthread -shared\
			   -Wl,--version-script,"$(MATLAB_ROOT)extern/lib/glnxa64/mexFunction.map" \
			   -Wl,--no-undefined -lblas -llapack
			   

			   
CXX_LIBS = -Wl,--no-undefined -Wl,-rpath-link,"$(MATLAB_ROOT)bin/glnxa64" \
		   -L$(MATLAB_ROOT)bin/glnxa64 -lmx -lmex -lmat -lm -fopenmp 
		   
		   

########################## Mesh ###########################		   
domWrapper    = src
domWrapperOut = class/tracer/private

$(domWrapper)/dom.o: $(domWrapper)/Tracer.cpp $(domWrapper)/Tracer.h
	$(CXX) -c $(CXX_INCLUDE) $(CXX_FLAGS) $< -o $@

$(domWrapperOut)/TracerWrapper.mexa64: $(domWrapper)/dom.o
	$(CXX) $(MATLAB_LINKS) -o $@ $< $(CXX_LIBS) && rm $(domWrapper)/dom.o
###########################################################

all:$(domWrapperOut)/TracerWrapper.mexa64

clean:
	rm -f $(domWrapperOut)/TracerWrapper.mexa64 