0;136;0cCC=gcc
CC+= -g
CFLAGS=-c -Wall
LDFLAGS=
SOURCES=Fit.C
OBJECTS=$(SOURCES:.C=.o)
EXECUTABLE=Fit

CFLAGS += $(shell $(ROOTSYS)/bin/root-config --cflags)
LDFLAGS += $(shell $(ROOTSYS)/bin/root-config --libs) -lMinuit -lMinuit2
LDFLAGS += -lTreePlayer -lstdc++

all: $(SOURCES) $(EXECUTABLE)

depend: .depend

.depend: $(SOURCES)
	rm -f ./.depend
	$(CC) $(CFLAGS) -MM $^ >> ./.depend;

include .depend

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

.C.o:
	$(CC) $(CFLAGS) $< -o $@

clean:
	rm -f *.o; rm -f */*.o; rm -f .depend; rm $(EXECUTABLE); rm -f out*.dat
