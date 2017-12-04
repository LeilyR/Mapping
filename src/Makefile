

CXX = g++ -g -std=c++0x
CFLAGS = -Wall -fmessage-length=50 -fopenmp #-march=bdver2
CXXFLAGS = -Wall -fmessage-length=50 -fopenmp -pedantic #-march=bdver2
#OPTFLAGS = -O3 -DNDEBUG
#OPTFLAGS = -O3 -DNDEBUG 
DEBUGFLAGS = -g3 -O0 -gdwarf-2 

DEBUG = 1
#OPT = 1

ifdef DEBUG
CXXFLAGS += $(DEBUGFLAGS)
endif

ifdef OPT
CXXFLAGS += $(OPTFLAGS)
endif

CXXFLAGS += -Wlong-long


#STATGEN_INC = libStatGen-1.0.13/include
#INCLUDEARGS = -I$(STATGEN_INC)
#LIBARGS =   apclusterunix64.so -L../dlib_build -ldlib  libStatGen-1.0.13/libStatGen.a -lz


# include libStatGen-1.0.13/Makefiles/Makefile.include

PROGS =  mapping

OBJS = data.o needleman_wunsch.o map_read.o dijkstra.o ref_graph.o test_mapping.o simple_NW_algo.o fibonacci_heap.o

all: $(PROGS)




mapping: main.cpp $(OBJS) pw_alignment.hpp data.hpp needleman_wunsch.hpp simple_NW_algo.hpp map_read.hpp ref_graph.hpp test_mapping.hpp alignment_index.hpp alignment_index.cpp intervals.hpp intervals.cpp needleman_wunsch.cpp dynamic_mc.hpp dynamic_mc.cpp fibonacci_heap.hpp fibonacci_heap.cpp
	$(CXX) -o mapping $(CXXFLAGS) $(INCLUDEARGS) main.cpp $(OBJS) $(LIBARGS)
	
data.o: data.cpp pw_alignment.hpp data.hpp
	$(CXX) -o data.o $(CXXFLAGS) $(INCLUDEARGS) -c data.cpp 


#dynamic_mc.o: dynamic_mc.cpp dynamic_mc.hpp  data.hpp
#	$(CXX) -o dynamic_mc.o $(CXXFLAGS) $(INCLUDEARGS) -c dynamic_mc.cpp

needleman_wunsch.o: needleman_wunsch.cpp needleman_wunsch.hpp
	$(CXX) -o needleman_wunsch.o $(CXXFLAGS) $(INCLUDEARGS) -c needleman_wunsch.cpp
dijkstra.o: dijkstra.hpp dijkstra.cpp
	$(CXX) -o dijkstra.o $(CXXFLAGS) $(INCLUDEARGS) -c dijkstra.cpp


map_read.o: map_read.cpp map_read.hpp 
	$(CXX) -o map_read.o $(CXXFLAGS) $(INCLUDEARGS) -c map_read.cpp

test_mapping.o: test_mapping.cpp test_mapping.hpp
	$(CXX) -o test_mapping.o $(CXXFLAGS) $(INCLUDEARGS) -c test_mapping.cpp


ref_graph.o: ref_graph.cpp ref_graph.hpp
	$(CXX) -o ref_graph.o $(CXXFLAGS) $(INCLUDEARGS) -c ref_graph.cpp

simple_NW_algo.o: simple_NW_algo.cpp simple_NW_algo.hpp
	$(CXX) -o simple_NW_algo.o $(CXXFLAGS) $(INCLUDEARGS) -c simple_NW_algo.cpp

fibonacci_heap.o: fibonacci_heap.cpp fibonacci_heap.hpp
	$(CXX) -o fibonacci_heap.o $(CXXFLAGS) $(INCLUDEARGS) -c fibonacci_heap.cpp

clean:
	@-rm -f $(PROGS) *.o  *.a *~





