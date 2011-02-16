PUB	= /home/linux/ieng6/cs160w/public
include $(PUB)/Arch/arch.gnu.omp
# include $(PUB)/Arch/arch.gnu.pthread

WARNINGS += -Wall -pedantic

# OPTIMIZATION = -pg -g

# Uncomment this line if your local installation of pthreads
# does not support Barriers 

# If your local installation of pthreads
# does not support Barriers 
# specify pbarrier=1 on the "make" command line
ifeq ($(pbarrier), 1)
    C++FLAGS += -DPTHREAD_BARRIER
    CFLAGS   += -DPTHREAD_BARRIER
else
endif


# This generates output about how the
# compiler vectorized the code
# We  suggest using level 2 (the integer after "verbose=")
# See the gcc manual for the other levels of output: levels 1-7
# http://gcc.gnu.org/onlinedocs/gcc-4.4.2/gcc/Debugging-Options.html#Debugging-Options
OPTIMIZATION += -ftree-vectorizer-verbose=2
# OPTIMIZATION += -ftree-vectorize


# Add symbol table information for gdb
ifeq ($(debug), 1)
	CFLAGS += -g
	C++FLAGS += -g
endif


app:	lu


OBJ	= lu.o elim.o verify.o init.o cmdLine.o Report.o Timer.o

lu:	$(OBJ)
	$(C++LINK) $(LDFLAGS) -o $@ $(OBJ)  $(LDLIBS)

clean:	
	$(RM) *.o lu
	$(RM) core
