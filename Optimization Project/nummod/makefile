_processor = $(shell uname -p)

MOL_VERSION = 0.0.4
MOL_INCLUDE ="\"mol.$(MOL_VERSION).h\""
CPPFLAGS := -D _MOL_VERSION_="\"$(MOL_VERSION)\"" -D _MOL_INCLUDE_=$(MOL_INCLUDE) $(CPPFLAGS)
CPPFLAGS := -D ATOM_PRM="\"atom.$(MOL_VERSION).prm\"" -static $(CPPFLAGS) 
CC=gcc

LDFLAGS = -L../lib
CPPFLAGS := -I../include $(CPPFLAGS)



CFLAGS = -O3 -Wall -W -Wshadow -Wpointer-arith -Wcast-qual
LIBS := -lmol.$(MOL_VERSION)  -lm $(LIBS)

OBJS = nummod.o

nummod:	$(OBJS)
	$(CC) $(CFLAGS) $(CPPFLAGS) $(LDFLAGS) $(OBJS) $(LIBS) -o nummod

%.o: %.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -c $<

all:
	nummod

clean:
	$(RM) $(OBJS)
	$(RM) nummod
