#
# 'make depend' uses makedepend to automatically generate dependencies
# (dependencies are added to end of Makefile)
# 'make' build executable file 'mycc'
# 'make clean' removes all .o and executable files
#
# define the C compiler to use
CC = g++

# define any compile-time flags
CFLAGS = -Wall -DHAVE_INLINE

HOSTNAME = $(shell if [ -d /home/shyamg/local/lib ]; then echo spud; else echo nono; fi)
# define any directories containing header files other than /usr/include
#
INCLUDES =
LFLAGS = 
ifeq ($(findstring spud, ${HOSTNAME}), spud)
	INCLUDES += -I/home/shyamg/local/include
	LFLAGS += -L/home/shyamg/local/lib
endif

# define library paths in addition to /usr/lib
# if I wanted to include libraries not in /usr/lib I'd specify
# their path using -Lpath, something like:

# define any libraries to link into executable:
# if I want to link in libraries (libx.so or libx.a) I use the -llibname
# option, something like (this will link in libmylib.so and libm.so:
LIBS = -lgsl -lgslcblas -lnlopt -lm

# define the CPP source files
SRCS = migrate.cpp main.cpp

# define the C object files
#
# This uses Suffix Replacement within a macro:
# $(name:string1=string2)
# For each word in 'name' replace 'string1' with 'string2'
# Below we are replacing the suffix .c of all words in the macro SRCS
# with the .o suffix
#
OBJS = $(SRCS:.cpp=.o)

# define the executable file
MAIN = migrate

#
# The following part of the makefile is generic; it can be used to
# build any executable just by changing the definitions above and by
# deleting dependencies appended to the file from 'make depend'
#

.PHONY: depend clean debug debuglocal local all

all: $(MAIN)
	@echo Migrate program has been compiled

$(MAIN): $(OBJS)
	$(CC) $(CFLAGS) $(INCLUDES) $(LFLAGS) -static -o $(MAIN) $(OBJS) $(LIBS)

# this is a suffix replacement rule for building .o's from .c's
# it uses automatic variables $<: the name of the prerequisite of
# the rule(a .c file) and $@: the name of the target of the rule (a .o file)
# (see the gnu make manual section about automatic variables)
.cpp.o:
	$(CC) $(CFLAGS) $(INCLUDES) -c $< -o $@

debug: CFLAGS += -g -DDEBUG
debug: $(MAIN)

debuglocal: CFLAGS += -g -DDEBUG -DLOCAL
debuglocal: $(MAIN)

local: CFLAGS += -DLOCAL
local: $(MAIN)

test: CFLAGS += -DDEBUG -DTEST
test: $(MAIN)

clean:
	rm -f *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^
# DO NOT DELETE

migrate.o: migrate.h
main.o: migrate.h
