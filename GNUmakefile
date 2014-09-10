# ===========================================================================
#  Makefile superMC                                    Chun Shen Mar. 19, 2013
# ===========================================================================
##
##  Environments :	MAIN	= 	main sourcefile	
##
##  Usage : 	(g)make	[all]		compile the whole project		
##			install		make all and copy binary to $INSTPATH
##			clean		remove objectfiles in obj_$TYPE 
##			distclean	remove all objectsfiles and binaries
##  

CC := g++
CFLAGS := -g $(shell gsl-config --cflags)

RM		=	rm -f
O               =       .o
LDFLAGS         =       $(CFLAGS) $(shell gsl-config --libs)
SYSTEMFILES     =       $(SRCGNU)

# --------------- Files involved ------------------

ifeq "$(MAIN)" ""
MAIN		=	nucleusParametersTuner.e
endif

SRC		=	main.cpp Nucleus.cpp ParameterReader.cpp

INC		= 	Nucleus.h ParameterReader.h Particle.h


# -------------------------------------------------

OBJDIR		=	obj
SRCFILES 	= 	$(SRC) $(INC) GNUmakefile
OBJECTS		=	$(addprefix $(OBJDIR)/, $(addsuffix $O, \
			$(basename $(SRC))))
TARGET		=	$(MAIN)
INSTPATH	=	$(HOME)/local/bin

# --------------- Pattern rules -------------------

$(OBJDIR)/%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

%.cpp:
	if [ -f $@ ] ; then touch $@ ; else false ; fi

# -------------------------------------------------

.PHONY:		all mkobjdir clean distclean install

all:		mkobjdir $(TARGET)

help:
		@grep '^##' GNUmakefile

mkobjdir:	
		-@mkdir -p $(OBJDIR)

$(TARGET):	$(OBJECTS)	
		$(CC) $(OBJECTS) -o $(TARGET) $(LDFLAGS) 
#		strip $(TARGET)

clean:		
		-rm $(OBJECTS)

distclean:	
		-rm $(TARGET)
		-rm -r obj

install:	$(TARGET)
		cp $(TARGET) $(INSTPATH)

# --------------- Dependencies -------------------
./main.cpp: Nucleus.h ParameterReader.h
./Nucleus.cpp: Nucleus.h ParameterReader.h
./ParameterReader.cpp: ParameterReader.h
