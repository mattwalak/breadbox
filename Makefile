all: main.cpp bread

SOURCES    = main.cpp
OBJECTS    = $(SOURCES:.cpp=.o)

.cpp.o:
	g++ -w -O3 -c $< -o $@

bread: main.o
	g++ $(OBJECTS) $(LDFLAGS) -o $@

clean:
	rm -f *.o

