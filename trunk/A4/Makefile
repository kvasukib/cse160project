PUB	= /home/linux/ieng6/cs160w/public/
include $(PUB)/Arch/arch.gnu.mpi

WARNINGS += -Wall -pedantic

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


# Add symbol table information for gdb/cachegrind
ifeq ($(debug), 1)
	CFLAGS += -g
	C++FLAGS += -g
endif



# Add symbol table information for gprof
ifeq ($(gprof), 1)
	CFLAGS += -g -pg
	C++FLAGS += -g -pg
	C++FLAGS += -g -pg
	LDFLAGS += -g -pg
endif

app:		apf

OBJECTS = apf.o solve.o splot.o cmdLine.o Report.o printMat.o 

apf:	        $(OBJECTS) 
		$(C++LINK) $(LDFLAGS) -o $@ $(OBJECTS)  $(LDLIBS)

clean:	
	$(RM) *.o apf;
	$(RM) core;
