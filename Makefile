CXX     = g++
CC     = gcc
#CFLAGS = -Wall -Wextra -O2 -fopenmp
CFLAGS = -Wall -Wextra -Og -fopenmp -ggdb
LDLIBS = -lm 

TGTS = heated-plate-parallel pi sequence_alignment

all: $(TGTS)

.PHONY: clean
clean:
	rm -f $(TGTS) *.o


heated-plate-parallel: heated-plate-parallel.o
	$(CXX) $(CFLAGS) -o heated-plate-parallel heated-plate-parallel.o $(LDLIBS)

heated-plate-parallel.o: heated-plate-parallel.cpp
	$(CXX) $(CFLAGS) -c heated-plate-parallel.cpp

pi: pi.o
	$(CC) $(CFLAGS) -o pi pi.o $(LDLIBS)

pi.o: pi.c
	$(CC) $(CFLAGS) -c pi.c

sequence_alignment: sequence_alignment.o
	$(CC) $(CFLAGS) -o sequence_alignment sequence_alignment.o $(LDLIBS)

sequence_alignment.o: sequence_alignment.c
	$(CC) $(CFLAGS) -c sequence_alignment.c

