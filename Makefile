# SRC_DIR := trident_smite
# OBJ_DIR := 
# SRC_FILES := $(wildcard $(SRC_DIR)/*.cxx)
# OBJ_FILES := $(patsubst $(SRC_DIR)/%.cxx,$(OBJ_DIR)/%.o,$(SRC_FILES))
# LDFLAGS := -lm -lcuba -I /usr/include/ ./$(SRC_DIR) -L /usr/lib/ -lboost_program_options
# CPPFLAGS := gen_SM.cxx integrator.cxx integrands.cxx cross_sections.cxx observables.cxx constants.cxx
# CXXFLAGS := -std=c++11 


# # $source/cross_sections.cxx $source/integrator.cxx $source/integrands.cxx $source/observables.cxx $source/integrator.cxx $source/constants.h
# make: $(OBJ_FILES)
# 	g++ $(LDFLAGS) -o $@ $^ 

# $(OBJ_DIR)/%.o: $(SRC_DIR)/%.cpp
# 	g++ $(CPPFLAGS) $(CXXFLAGS) -c -o $@ $<


#Compiler and Linker
CC          := g++

#The Target Binary Program
TARGET      := gen_SM

# Project folder
PROJ        := smite/

#The Directories, Source, Includes, Objects, Binary and Resources
SRCDIR      := $(PROJ)src
INCDIR      := $(PROJ)inc
BUILDDIR    := $(PROJ)obj
TARGETDIR   := ./
RESDIR      := $(PROJ)res
SRCEXT      := cxx
DEPEXT      := d
OBJEXT      := o

#Flags, Libraries and Includes
CFLAGS      := -O3 -g --std=c++17 #-Wall 
LIB         := -lm -lcuba -lz -L /usr/lib/ -lboost_program_options
INC         := -I $(INCDIR) -I /usr/include 
INCDEP      := -I $(INCDIR)

#---------------------------------------------------------------------------------
#DO NOT EDIT BELOW THIS LINE
#---------------------------------------------------------------------------------
SOURCES     := $(shell find $(SRCDIR) -type f -name *.$(SRCEXT))
OBJECTS     := $(patsubst $(SRCDIR)/%,$(BUILDDIR)/%,$(SOURCES:.$(SRCEXT)=.$(OBJEXT)))

#Defauilt Make
all: resources $(TARGET)

#Remake
remake: cleaner all

#Copy Resources from Resources Directory to Target Directory
resources: directories
	@cp -r $(RESDIR)/* $(TARGETDIR)/

#Make the Directories
directories:
	@mkdir -p $(TARGETDIR)
	@mkdir -p $(BUILDDIR)

#Clean only Objecst
clean:
	@$(RM) -rf $(BUILDDIR)

#Full Clean, Objects and Binaries
cleaner: clean
	@$(RM) -rf $(TARGETDIR)

#Pull in dependency info for *existing* .o files
-include $(OBJECTS:.$(OBJEXT)=.$(DEPEXT))

#Link
$(TARGET): $(OBJECTS)
	$(CC) -o $(TARGETDIR)/$(TARGET) $^ $(LIB)

#Compile
$(BUILDDIR)/%.$(OBJEXT): $(SRCDIR)/%.$(SRCEXT)
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS) $(INC) -c -o $@ $<
	@$(CC) $(CFLAGS) $(INCDEP) -MM $(SRCDIR)/$*.$(SRCEXT) > $(BUILDDIR)/$*.$(DEPEXT)
	@cp -f $(BUILDDIR)/$*.$(DEPEXT) $(BUILDDIR)/$*.$(DEPEXT).tmp
	@sed -e 's|.*:|$(BUILDDIR)/$*.$(OBJEXT):|' < $(BUILDDIR)/$*.$(DEPEXT).tmp > $(BUILDDIR)/$*.$(DEPEXT)
	@sed -e 's/.*://' -e 's/\\$$//' < $(BUILDDIR)/$*.$(DEPEXT).tmp | fmt -1 | sed -e 's/^ *//' -e 's/$$/:/' >> $(BUILDDIR)/$*.$(DEPEXT)
	@rm -f $(BUILDDIR)/$*.$(DEPEXT).tmp

#Non-File Targets
.PHONY: all remake clean cleaner resources