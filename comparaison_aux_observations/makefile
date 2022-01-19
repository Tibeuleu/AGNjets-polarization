CC=g++
CFLAGS=-c -Wall --pedantic
LDFLAGS=-Wall --pedantic
SOURCES=functions.cpp comparaison.cpp
OBJECTS=functions.o comparaison.o
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
comparaison.o: functions.h
