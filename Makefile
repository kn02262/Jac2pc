FLAGS = -O4 -march=native -ftree-vectorize -funroll-loops -std=c++11 -Wall

LIBS = /usr/local/lib/libntl.a -lgmp -pthread

CXX = g++

objects = elemQuo.o resultantsumroots.o ZZ_pEJac2F.o ZZ_pXResultant.o elltorsion.o schoof.o ZZ_pJac2.o fastinterp.o vec_pair_long_long.o ZZ_pJac2F.o ZZ_pEJac2.o ZZ_pXCantorpoly.o cpuperf.o RW2dim.o LMPMCT_multithread.o LMPMCT_onethread.o

all: $(objects) main LMPMCT

$(objects): %.o: %.c %.h
	$(CXX) -c $(FLAGS) $< -o $@

main: main.c $(objects)
	$(CXX) -o main main.c $(objects) $(FLAGS) $(LIBS) 

RandCurveGen: Random_curve_generator.c $(objects)
	$(CXX) -o Random_curve_generator Random_curve_generator.c $(objects) $(FLAGS) $(LIBS) 

LMPMCT: LMPMCT.c LMPMCT_multithread.c LMPMCT_onethread.c $(objects)
	$(CXX) -o LMPMCT LMPMCT.c $(objects) $(FLAGS) $(LIBS)
	
clean:
	rm -f *.o main LMPMCT testLauter
