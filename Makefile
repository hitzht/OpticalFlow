SHELL='bash'
#
# Bridges - PSC
#
# Intel Compilers are loaded by default
# You will need to specifically switch to GNU Modules
# With with `modules.sh`
#

CC = g++
MPCC = mpic++
OPENMP = -fopenmp
CFLAGS = -O3
LIBS =

cuCC = nvcc
cuCFLAGS = -std=c++11 -O3 -arch=compute_30 -code=sm_30
cuNVCCFLAGS = -std=c++11 -O3 -arch=compute_30 -code=sm_30

TARGETS = Matching_Ver1 Matching_Ver2 Matching_Ver3 Matching_Ver4 Matching_Ver5 Matching_Ver6 Matching_Ver7

all:	$(TARGETS)

Matching_Ver1: Matching_Ver1.cpp
	$(CC) -o $@ $(LIBS) Matching_Ver1.cpp
Matching_Ver2: Matching_Ver2.cpp
	$(CC) -o $@ $(LIBS) Matching_Ver2.cpp
Matching_Ver3: Matching_Ver3.cpp
	$(CC) -o $@ $(LIBS) Matching_Ver3.cpp
Matching_Ver4: Matching_Ver4.cu 
	$(cuCC) -o $@ $(LIBS) $(cuNVCCFLAGS) Matching_Ver4.cu
Matching_Ver5: Matching_Ver5.cpp
	$(CC) -o $@ $(LIBS) -lpthread Matching_Ver5.cpp
Matching_Ver6: Matching_Ver6.cpp
	$(CC) -o $@ $(LIBS) $(OPENMP) Matching_Ver6.cpp
Matching_Ver7: Matching_Ver7.cpp
	$(MPCC) -o $@ $(LIBS) $(MPILIBS) Matching_Ver7.cpp

clean:
	rm -f *.o $(TARGETS) *.stdout *.txt