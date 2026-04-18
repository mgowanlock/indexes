COMMONINCLUDES = params.h compare.h kdtree.h structs.h


###############################
#Host variables
CC=g++
LDFLAGS := -fopenmp -std=c++11 -O3 
###############################

all: main


main.o: main.cpp $(COMMONINCLUDES)
	$(CC) -c main.cpp $(LDFLAGS) -o main.o

main: main.o 
	$(CC) main.o -o main $(LDFLAGS)


clean:
	rm *.o
	rm main 
	


