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
#CXXFLAGS= -g
CXXFLAGS= 
BOOST_ROOT= /usr/local/include/
KONA_ROOT= /home/jehicken/Kona

# linker options
# NOTE: -L is for linking, -Wl,-rpath is for loading
#LDFLAGS= -lstdc++ -L$(KONA_ROOT) -Wl,-rpath,$(KONA_ROOT) -lkona
LDFLAGS= -lstdc++ -L$(BOOST_ROOT)/stage/lib -lboost_program_options \
	-L$(KONA_ROOT) -lkona

# options that DO NOT vary
ALL_CXXFLAGS= -I. $(CXXFLAGS) -I $(BOOST_ROOT) -I $(KONA_ROOT)/src

# source directory
SRC=../

# source, object, and executable file names
HEADERS= $(wildcard $(SRC)*.hpp)
SOURCES= $(wildcard $(SRC)*.cpp)
OBJS= $(SOURCES:.cpp=.o)
BINARIES= design.bin

# implicit rule
%.o : %.cpp $(HEADERS) Makefile
	@echo "Compiling \""$@"\" from \""$<"\""
	@$(CXX) $(CPPFLAGS) $(ALL_CXXFLAGS) -o $@ -c $<

default: all

all: $(BINARIES)

design.bin: $(OBJS) design.o Makefile
	@echo "Compiling \""$@"\" from \""$(OBJS)"\""
	@$(CXX) -o $@ design.o $(OBJS) $(LDFLAGS)

clean:
	@echo "deleting temporary, object, and binary files"
	@rm -f *~
	@rm -f $(BINARIES)
	@rm -f $(OBJS) *.o