#------------------------------------------------------------------------------
#                  Central Makefile for the program KERNELS
#------------------------------------------------------------------------------

# Mark Hoggard - September 2016.
# Run on Jupiter

# Source code originally from Sia Ghelichkhan
# NB. I couldn't get GFORTRAN to work. Instead using IFORT on Jupiter (Keith's machine) and building statically (all code included in executable)

# Set directory information
TARGETDIR = .
SRCDIR    = CODE/SRC
OBJDIR    = CODE/OBJ

# Set name executable
BIN = $(TARGETDIR)/KERNELS

# Compiler flags for IFORT (NB. Different if want to use GFORTRAN)
LK	= gfortran
FC	= gfortran
FFLAGS   = -fdefault-real-8 -fdefault-double-8 -fdefault-integer-8 -ffast-math
MOREFFLAGS =  -funroll-loops -march=native -cpp -Ofast -w
LFLAGS   = -llapack 

# List of object files - generate them if they don't exist using the builder
OBJ = $(OBJDIR)/driver.o $(OBJDIR)/calc.o $(OBJDIR)/matrix.o $(OBJDIR)/util.o

# Generate object directory, object files and the final executable
$(BIN) : $(OBJDIR) $(OBJ)
	 $(LK) -o $(BIN) $(OBJ) $(LFLAGS)

# Build objects
$(OBJDIR)/%.o : $(SRCDIR)/%.f $(SRCDIR)/grid.h $(SRCDIR)/parameter.h $(SRCDIR)/size.h
	        $(FC) -o $@ -c $< $(FFLAGS) $(MOREFFLAGS)

# Make the object directory if it doesn't exist
$(OBJDIR) : ;
	mkdir $(OBJDIR)

# Remove all object executables and the final executable
clean:
	$(RM) $(OBJDIR)/*.o
	$(RM) $(BIN)
	$(RM) -r $(OBJDIR)

