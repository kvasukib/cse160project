PUB	= /home/linux/ieng6/cs160w/public/
include $(PUB)/Arch/arch.gnu.pthread

WARNINGS += -Wall -pedantic

# Uncomment this line if your local installation of pthreads
# does not support Barriers 
# C++FLAGS += -DPTHREAD_BARRIER

# Serial Version (Provided)
OBJ	= Primes.o isPrime.o init.o Report.o Timer.o 

primes:	        $(OBJ)
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJ)  $(LDLIBS)


# Multithreaded Version
TOBJ	= thPrimes.o prime_thr.o init_thr.o Report.o isPrime.o parTimer.o Timer.o

th-primes:	$(TOBJ)
		$(C++LINK) $(LDFLAGS) -o $@ $(TOBJ)  $(LDLIBS)

clean:	
	$(RM) *.o primes th-primes
	$(RM) core
