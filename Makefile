# Define variable CC to be the compiler we want to use
CC = g++



TARGETS = clean AdaptiveController

# If no arguments are passed to make, it will attempt the 'c-example' target
default: AdaptiveController

# This runs the 'clean' and 'c-example' targets when 'make all' is run
all: $(TARGETS)

# This will construct the binary 'c-example'
# $^ = names of all the dependent files, deduped and with spaces
# $@ = complete name of the target
AdaptiveController: AdaptiveControl.cpp
	$(CC) $(CFLAGS) $^ `wx-config --cxxflags --libs std` -o $@

# $(RM) is the platform agnostic way to delete a file (here rm -f)
clean:
	$(RM) -rf AdaptiveController *~ *dSYM
