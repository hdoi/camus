# Makefile
F90= gfortran
LIB = 
OPT = -O3 -fopenmp

TARGET = camus
OBJS =  v5.o sub.o rand_machine.o
SRC =   sub.f90
MOD = rand_machine.f90 v5.f90
OBJ = $(addprefix obj/,$(SRC:.f90=.o))
OBJM= $(addprefix obj/,$(MOD:.f90=.o))
OBJS= $(OBJM) $(OBJ)
OBJDIR= obj
.SUFFIXES:
.SUFFIXES: .f90 .o

$(TARGET) : $(OBJS)
	$(F90) -o bin/$(TARGET) -I $(OBJDIR) $(OPT) $(OBJS) $(LIB)

$(OBJDIR)/%.o : src/%.f90
	$(F90) $(OPT) -I $(OBJDIR) -J $(OBJDIR) -o $@ -c $<

$(OBJ): $(addprefix src/,$(MOD))

clean :
	rm -f $(OBJDIR)/*.mod $(OBJDIR)/*.o bin/$(TARGET)

