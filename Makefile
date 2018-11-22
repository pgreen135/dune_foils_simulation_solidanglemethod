CXXFLAGS=-std=c++11 $(shell root-config --cflags)
LIBS=$(shell root-config --libs) -lMathMore

run : analyze_light_histo
			@echo "Finished Compiling..."
			@echo "To run: ./analyze_light_histo"

analyze_light_histo : analyze_light_histo.o utility_functions.o timeparam.o solid_angle_functions.o

	g++ -o $@ $^ ${LIBS}

%.o : %.cc
	g++ ${CXXFLAGS} -o $@ -c $^
