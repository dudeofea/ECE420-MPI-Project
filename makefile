FLAGS=--std=c99 -Wall

all:
	$(CC) $(FLAGS) datatrim.c -o datatrim
	$(CC) $(FLAGS) serialcalc.c -o serialcalc -lm
	mpicc $(FLAGS) mpicalc.c -o mpicalc -lm
run:
	mpirun -np 8 ./mpicalc