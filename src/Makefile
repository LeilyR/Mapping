

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


STATGEN_INC = libStatGen-1.0.13/include
INCLUDEARGS = -I$(STATGEN_INC)
LIBARGS =   apclusterunix64.so -L../dlib_build -ldlib  libStatGen-1.0.13/libStatGen.a -lz


# include libStatGen-1.0.13/Makefiles/Makefile.include

PROGS =  mapping

OBJS = data.o needleman_wunsch.o map_read.o dijkstra.o 

all: $(PROGS)




mapping: main.cpp $(OBJS) pw_alignment.hpp data.hpp needleman_wunsch.hpp map_read.hpp alignment_index.hpp alignment_index.cpp intervals.hpp intervals.cpp needleman_wunsch.cpp dynamic_mc.hpp dynamic_mc.cpp
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

#connected_component.o: connected_component.cpp connected_component.hpp
#	$(CXX) -o connected_component.o $(CXXFLAGS) $(INCLUDEARGS) -c connected_component.cpp


#model.o: model.cpp model.hpp
#        $(CXX) -o model.o $(CXXFLAGS) $(INCLUDEARGS) -c model.cpp

clean:
	@-rm -f $(PROGS) *.o  *.a *~





