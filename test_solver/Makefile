# makefile for quasi_1d_euler

SHELL = /bin/sh

.SUFFIXES:
.SUFFIXES: .cpp .o
.PHONY: default clean all solver

# compiler
CXX= gcc 

# c preprocessor options
CPPFLAGS= -cpp -DNDEBUG -DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS -DBOOST_UBLAS_TYPE_CHECK=0
#CPPFLAGS= -cpp -DDEBUG -DBOOST_UBLAS_ENABLE_PROXY_SHORTCUTS

# compiler options that may vary (user can change)
CXXFLAGS= 
BOOST_ROOT= $(BOOST_HOME)
KONA_ROOT= $(KONA_HOME)

# linker options
# NOTE: -L is for linking, -Wl,-rpath is for loading
#LDFLAGS= -lstdc++ -L$(KONA_ROOT) -Wl,-rpath,$(KONA_ROOT) -lkona
LDFLAGS= -lm -lstdc++ -L$(BOOST_ROOT)/stage/lib -lboost_program_options \
	-L$(KONA_ROOT) -lkona

# options that DO NOT vary
ALL_CXXFLAGS= -I. $(CXXFLAGS) -I $(BOOST_ROOT) -I $(KONA_ROOT)/src

# source directory
SRC=../

# source, object, and executable file names
HEADERS= $(wildcard $(SRC)*.hpp)
SOURCES= $(wildcard $(SRC)*.cpp)
OBJS= $(SOURCES:.cpp=.o)
BINARIES= verify.bin solver.bin

# implicit rule
%.o : %.cpp $(HEADERS) Makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CPPFLAGS) $(ALL_CXXFLAGS) -o $@ -c $<

default: all

all: $(BINARIES)

verify.bin: $(OBJS) verify_accuracy.o Makefile
	@echo "Compiling \""$@"\" from \""$(OBJS)"\""
	@$(CXX) -o $@ verify_accuracy.o $(OBJS) $(LDFLAGS)

solver.bin: $(OBJS) test_quasi_1d_euler.o Makefile
	@echo "Compiling \""$@"\" from \""$(OBJS)"\""
	@$(CXX) -o $@ test_quasi_1d_euler.o $(OBJS) $(LDFLAGS) 

clean:
	@echo "deleting temporary, object, and binary files"
	@rm -f *~
	@rm -f $(BINARIES)
	@rm -f $(OBJS) *.o