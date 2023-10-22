# Set the C compiler
#CC            =

HDRS	        = dg.h dg1d.h

LD            = $(CC)

MAKEFILE      = makefile

OBJS          = main.o init.o shape.o gauss.o flux.o atgauss.o \
                save.o update.o div.o result.o project.o time.o

PROGRAM       = dg

all:            $(PROGRAM)

$(PROGRAM):     $(OBJS)
				        $(LD) -o $(PROGRAM) $(OBJS) -lm

clean:;         rm -f $(OBJS) core $(PROGRAM)
