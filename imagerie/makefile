CC=g++
CFLAGS=-c -Wall --pedantic
LDFLAGS=-Wall --pedantic
SOURCES=synchrotron.cpp functions.cpp simulation_image.cpp
OBJECTS=synchrotron.o functions.o simulation_image.o
EXECUTABLE=projet

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ ; ./$(EXECUTABLE) ; make clean ;

clean:
	$(RM) *.o $(EXECUTABLE)

depend:
	makedepend -Y $(SOURCES)
# DO NOT DELETE

synchrotron.o: functions.h
functions.o: functions.h boost/math/special_functions/bessel.hpp
simulation_image.o: functions.h synchrotron.h
