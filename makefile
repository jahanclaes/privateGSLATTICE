# Compiler
CC = mpiicpc
# Optimization level (-xhost only works for intel compilers)
CC_OPT = -O3 -xhost
# Compilation flags (Adding -std=c++11 will cause trouble for intel compiler)
CC_FLAGS = -DUSE_MPI -openmp -pthread -mkl=parallel
# Include directories
CC_INC = -I/media/bkcgroup/CoolStuff/Install/include/
# Linking directories
CC_LIB = -L/media/bkcgroup/CoolStuff/Install/lib/ -L/opt/local/lib/ -lsprng

# File names
EXEC = run
SOURCES = $(wildcard ./*.cc ./cPEPS_Base_Class/*.cc ./Communication/*.cc)
OBJECTS = $(SOURCES:.cc=.o)

# Main target
$(EXEC): $(OBJECTS)
	$(CC) $(OBJECTS) $(CC_OPT) $(CC_FLAGS) $(CC_INC) $(CC_LIB) -o $(EXEC)

# To obtain object files
%.o: %.cc
	$(CC) -c $(CC_OPT) $(CC_FLAGS) $(CC_INC) $(CC_LIB) $< -o $@

# To remove generated files
clean:
	rm -f *.o cPEPS_Base_Class/*.o PEPS_Base_Class/*.o Communication/*.o
