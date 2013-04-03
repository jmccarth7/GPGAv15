PROG =	GPCODE_test

SRCS =	0GPCODE_test.f90  \
	GP_Clean_Tree_Nodes.f90 GP_Elitists.f90 \
	GP_Fitness_Proportionate_Asexual_Reproduction.f90 GP_Mutations.f90 \
	GP_Tournament_Style_Sexual_Reproduction.f90 GP_Tree_Build.f90 \
	GP_Tree_Swap.f90 GPCODE_modules.f90  \
	random_real.f90

OBJS =	0GPCODE_test.o  \
	GP_Clean_Tree_Nodes.o GP_Elitists.o \
	GP_Fitness_Proportionate_Asexual_Reproduction.o GP_Mutations.o \
	GP_Tournament_Style_Sexual_Reproduction.o GP_Tree_Build.o \
	GP_Tree_Swap.o GPCODE_modules.o  \
	random_real.o

LIBS =	

CC = cc
CFLAGS = -O
FC = g95
FFLAGS = -g
F90 = g95
F90FLAGS = -g
#LDFLAGS = -lSystemStubs
LDFLAGS = -Wl,-no_pie
LIBS= -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.8.sdk/usr/lib -L/Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.7.sdk/usr/lib  -L/Developer/SDKs/MacOSX10.6.sdk/usr/lib



all: $(PROG)

$(PROG): $(OBJS)
	$(F90) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

new: $(PROG)

	rm -f  $(OBJS) *.mod

clean:
	rm -f $(PROG) $(OBJS) *.mod

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

0GPCODE_test.o: GPCODE_modules.o GPCODE_modules.f90 \
	Lotka_Volterra_Example_Set_Up.f901
GP_Clean_Tree_Nodes.o: GPCODE_modules.o
GP_Elitists.o: GPCODE_modules.o
GP_Fitness_Proportionate_Asexual_Reproduction.o: GPCODE_modules.o
GP_Mutations.o: GPCODE_modules.o
GP_Tournament_Style_Sexual_Reproduction.o: GPCODE_modules.o
GP_Tree_Build.o: GPCODE_modules.o
GP_Tree_Swap.o: GPCODE_modules.o
