CFLAGS = -std=c17 -Wall -O3 -fopenmp
CXXFLAGS = -std=c++17 -Wall -O3 -fopenmp
LDFLAGS = -lfmt -lpthread

GENERATOR_SOURCES = graph_generator.c make_graph.c splittable_mrg.c utils.c
GENERATOR_HEADERS = graph_generator.h make_graph.h mod_arith_32bit.h mod_arith_64bit.h mod_arith.h splittable_mrg.h utils.h mrg_transitions.c
GENERATOR_OBJECTS = $(GENERATOR_SOURCES:.c=.o)

all: generator_omp


$(GENERATOR_OBJECTS): $(GENERATOR_SOURCES) $(GENERATOR_HEADERS)
	gcc $(CFLAGS) -c $(GENERATOR_SOURCES)

generator_omp: generator_omp.cpp $(GENERATOR_OBJECTS)
	g++ $(CXXFLAGS) -o generator_omp generator_omp.cpp $(LDFLAGS) $(GENERATOR_OBJECTS)

clean:
	rm -f *.o
	rm -f generator_omp

.PHONY: all obj clean
