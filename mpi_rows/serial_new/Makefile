CC := gcc
CCFLAGS := -O3 -fopenmp -Wall -g
INCLUDE := .
#--- VPATH contains the search path for the src files !
# VPATH := .
#--- Now find & list all c files on directories
SOURCES := $(wildcard ./*.c)
#--- Object files shouldn’t have the full path of
#--- SOURCES when we replace the .cpp by .o !
OBJECTS := $(addprefix ./obj/, $(notdir $(SOURCES:.c=.o)))
#####################################################
all: ./foxes-rabbits
./foxes-rabbits: $(OBJECTS)
	@echo $(OBJECTS)
	$(CC) $(CCFLAGS) -o $@ $^ -I $(INCLUDE)
#--- Now every .c file can easily be found from VPATH !
./obj/%.o: %.c
	mkdir -p obj
	$(CC) $(CCFLAGS) -c $^ -o $@ -I $(INCLUDE)
clean:
	rm ./obj/*.o ./foxes-rabbits
