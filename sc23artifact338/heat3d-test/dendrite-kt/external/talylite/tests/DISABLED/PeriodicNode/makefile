include ../../makefile.config

INCLUDES = -I. $(TALYFEM_INCLUDES)
LIBS = $(TALYFEM_LIBS)

VPATH = src

####### Targets.

TARGET	= pbc_node_test

##### Build rules.
$(TARGET): pbc_node_test.o 
	$(CC) -O -o $@ pbc_node_test.o $(LIBS)

all: $(TARGET)
	rm *.o

clean:
	rm -f *.o $(TARGET) tests/out*.plt tests/dd*.plt tests/config.txt

.PHONY: test

test:
	@cd tests; ./run_test.py
  
##### Implicit rules.
.SUFFIXES: .cpp 

.cpp.o:
	$(CC) -c $(CC_FLAGS) $(TALYFEM_FLAGS) $(INCLUDES) $<
