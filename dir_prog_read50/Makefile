PROG =	prog_read50

SRCS =	prog_read50.f90

OBJS =	prog_read50.o

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

