MPICC = mpic++
EXEC = mc_pi 

all: ${EXEC}

mc_pi: monte_carlo_pi_parallel.cpp
	${MPICC} -o mc_pi monte_carlo_pi_parallel.cpp

clean:
	rm ${EXEC}
