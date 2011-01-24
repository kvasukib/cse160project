PUB	= /home/linux/ieng6/cs160w/public/
# PUB =    /Users/baden/lib
include $(PUB)/Arch/arch.gnu.pthread

WARNINGS += -Wall -pedantic

# Uncomment this line if your local installation of pthreads
# does not support Barriers 

# C++FLAGS += -DPTHREAD_BARRIER
# If you want to compile for single precision,
# uncomment the next line
# ARITH = -DFLOAT

# This generates output about how the
# compiler vectorized the code
# We  suggest using level 2 (the integer after "verbose=")
# See the gcc manual for the other levels of output: levels 1-7
# http://gcc.gnu.org/onlinedocs/gcc-4.4.2/gcc/Debugging-Options.html#Debugging-Options
#
OPTIMIZATION += -ftree-vectorizer-verbose=2

C++FLAGS += $(ARITH)
CFLAGS += $(ARITH)

#DEBUG += -DDEBUG

#Uncomment the next line  if you want to add symbol table
# information for gdb
#DEBUG += -g


app:		apf

OBJECTS = apf.o solve.o splot.o cmdLine.o Report.o printMat.o Timer.o

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o apf;
	$(RM) core;
