CC=g++
CFLAGS=-c -Wall --pedantic
LDFLAGS=-Wall --pedantic
SOURCES=functions.cpp emission_ordonnee.cpp
OBJECTS=functions.o emission_ordonnee.o
EXECUTABLE=projet

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(LDFLAGS) $^ -o $@ ; ./$(EXECUTABLE) ; make clean ;

clean:
	$(RM) *.o $(EXECUTABLE)

depend:
	makedepend -Y $(SOURCES)
# DO NOT DELETE

functions.o: functions.h
emission_ordonnee.o: functions.h
