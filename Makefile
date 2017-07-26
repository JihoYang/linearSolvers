#Include files
SOURCES=helper.c main.c linear_solvers.c

#Compiler
#--------
CC = gcc-7
CFLAGS = -std=gnu99 -fstrict-overflow -Werror -Wshadow -Wstrict-overflow=4 -pedantic

#Linker flags
#------------
LDFLAGS= -lm

OBJECTS=$(SOURCES:.c=.o)
EXECUTABLE=sim

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CC) $(OBJECTS) -o $@ $(LDFLAGS)

clean:
	rm -f $(OBJECTS) $(EXECUTABLE)

$(OBJECTS): %.o : %.c
	$(CC) $(CFLAGS) -c $< -o $@
