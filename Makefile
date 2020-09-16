# 1. Put this file in the same folder as your 'driver' code 
#    (the code containing the 'main' function).

# 2. Edit LIBRARY_DIR to point at the location of your ITensor Library
#    source folder (this is the folder that has options.mk in it).
#    Also, edit TDVP_DIR to point at the location of the TDVP source files.
#LIBRARY_DIR=$(HOME)/software/itensor
LIBRARY_DIR=/public1/home/sc51526/soft/itensor
TDVP_DIR=/public1/home/sc51526/soft/TDVP-master

# 3. If your 'main' function is in a file called 'myappname.cc', then
#    set APP to 'myappname'. Running 'make' will compile the app.
#    Running 'make debug' will make a program called 'myappname-g'
#    which includes debugging symbols and can be used in gdb (Gnu debugger);
APP=mbl_tdvp

# 4. Add any headers your program depends on here. The make program
#    will auto-detect if these headers have changed and recompile your app.
HEADERS=$(TDVP_DIR)/tdvp.h

# 5. For any additional .cc (source) files making up your project,
#    add their full filenames here.
#CCFILES=$(APP).cc

#################################################################
#################################################################
#################################################################
#################################################################

include $(LIBRARY_DIR)/this_dir.mk
include $(LIBRARY_DIR)/options.mk

TENSOR_HEADERS=$(LIBRARY_DIR)/itensor/core.h

CCFLAGS+=-I$(TDVP_DIR)
CCGFLAGS+=-I$(TDVP_DIR)

CCFLAGS+=-fopenmp
CCGFLAGS+=-fopenmp


#Mappings --------------
#OBJECTS=$(patsubst %.cc,%.o,$(wildcard *.cc), $(CCFILES))
OBJECTS=$(patsubst %.cc,%.o,$(wildcard *.cc))

#Rules ------------------

%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCFLAGS) -o $@ $<

.debug_objs/%.o: %.cc $(HEADERS) $(TENSOR_HEADERS)
	$(CCCOM) -c $(CCGFLAGS) -o $@ $<

#Targets -----------------

build: $(APP)
debug: $(APP)-g

$(APP): $(OBJECTS) $(ITENSOR_LIBS)
	$(CCCOM) -o $(APP) $(CCFLAGS) $(OBJECTS) $(LIBFLAGS) ${LDFLAGS}

$(APP)-g: mkdebugdir $(GOBJECTS) $(ITENSOR_GLIBS)
	$(CCCOM) $(CCGFLAGS) $(GOBJECTS) -o $(APP)-g $(LIBGFLAGS)

clean:
	rm -fr .debug_objs *.o $(APP) $(APP)-g

mkdebugdir:
	mkdir -p .debug_objs

