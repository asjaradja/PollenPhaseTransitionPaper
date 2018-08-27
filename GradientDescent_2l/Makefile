CC = /usr/bin/g++
CC_FLAGS = -o 
C_FLAGS = -I/home/asja/data/wigxjpf-1.5/inc -I/home/asja/data/eigen
LD_FLAGS = -L/home/asja/data/wigxjpf-1.5/lib
LD_LIBS = -lwigxjpf -lm

all: main

main: main.cpp
	$(CC) -std=gnu++0x -O3 main.cpp $(CC_FLAGS) main $(C_FLAGS) $(LD_FLAGS) $(LD_LIBS)

clean: 
	rm main

#make or make main should compile your program
#make clean should remove it

