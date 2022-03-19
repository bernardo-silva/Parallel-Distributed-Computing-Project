CC := gcc
CCFLAGS := -O3 -fopenmp -Wall
INCLUDE := ../src
#--- VPATH contains the search path for the src files !
VPATH := ../src:../main
#--- Now find & list all c files on directories
SOURCES := $(wildcard ../src/*.c ../main/*.c)
#--- Object files shouldn’t have the full path of
#--- SOURCES when we replace the .cpp by .o !
OBJECTS := $(addprefix ../bin/, $(notdir $(SOURCES:.c=.o)))
#####################################################
all: ../bin/foxes-rabbits
../bin/foxes-rabbits: $(OBJECTS)
	@echo $(OBJECTS)
	$(CC) $(CCFLAGS) -o $@ $^ -I $(INCLUDE)
#--- Now every .c file can easily be found from VPATH !
../bin/%.o: %.c
	$(CC) $(CCFLAGS) -c $^ -o $@ -I $(INCLUDE)
clean:
	rm ../bin/*.o ../bin/foxes-rabbits
