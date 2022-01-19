CC=g++
CFLAGS=-c -Wall --pedantic
LDFLAGS=-Wall --pedantic
SOURCES=functions.cpp synchrotron.cpp emission_pola.cpp
OBJECTS=functions.o synchrotron.o emission_pola.o
EXECUTABLE=projet

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ ; ./$(EXECUTABLE) ; make clean ;

clean:
	$(RM) *.o $(EXECUTABLE)

depend:
	makedepend -Y $(SOURCES)
# DO NOT DELETE

functions.o: functions.h boost/math/special_functions/bessel.hpp
synchrotron.o: functions.h
emission_pola.o: functions.h synchrotron.h
