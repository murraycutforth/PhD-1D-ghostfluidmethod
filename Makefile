#
# Makefile
# 1D_GFM
#
# Based on the makefile at http://hiltmon.com/
#

# Compiler
CC := g++

# Folders
SRCDIR := source
BUILDDIR := objectfiles

# Target
TARGET := 1D_GFM.exe

# Code Lists
SOURCES := $(shell find $(SRCDIR) -type f -name *.cpp)
OBJECTS := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.cpp=.o))
EXTOBJS := ../exact_riemann_solver/exact_RS_idealgas.o

# Folder Lists
INCDIRS := include
INCLIST := -I include
BUILDLIST := $(BUILDDIR)
EXTINCLIST := -I ../exact_riemann_solver/

# Shared Compiler Flags
OPLEVEL := 
CFLAGS := -Wall -c -g -std=c++11 $(OPLEVEL)
INC := $(INCLIST)
TESTERFLAGS :=

# Linking Step
$(TARGET): $(OBJECTS)
	@mkdir -p $(BUILDLIST)
	@echo "Linking..."
	@echo "  Linking $(TARGET)"; $(CC) $^ $(EXTOBJS) $(OPLEVEL) -o $(TARGET)
	@echo "Success."

# Compilation Step
$(BUILDDIR)/%.o: $(SRCDIR)/%.cpp
	@mkdir -p $(BUILDLIST)	
	@echo "Compiling $<..."; $(CC) $(CFLAGS) $(INC) $(EXTINCLIST) -o $@ $<

# Test Environment
tester:
	@echo "Compiling maintest..."; $(CC) $(CFLAGS) $(TESTERFLAGS) $(INC) -o test/maintest.o test/maintest.cpp
	@echo "Linking maintest..."; $(CC) $(filter-out $(BUILDDIR)/main.o,$(OBJECTS)) test/maintest.o $(OPLEVEL) $(TESTERFLAGS) -o test/maintest
	@echo "Success."

clean:
	@echo "Cleaning $(TARGET)..."; rm -r $(BUILDDIR) $(TARGET) test/*.o test/maintest


