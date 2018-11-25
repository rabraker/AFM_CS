# ROOT=/home/arnold/sandbox/c_sandbox/plplot_sandbox

CC  = gcc
CFLAGS=
#CFLAGS = -Wall -I/usr/include -lgsl -lgslcblas -lm
#include makedefs
APP_SRC= test_data.c main.c SOS.c fit_ZPK.c plotting.c
# plplot_simple.c
APP_NAME = bode_fit

INCLUDE  = -I/usr/include/plplot/ -Icmpfit/ -Iinclude
LIB=-L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/plplot5.10.0 -Lcmpfit/

#
#
# addprefix: $(addprefix src/,foo bar) --> 'src/foo', 'src/bar'
# https://www.gnu.org/software/make/manual/make.html#index-addprefix

# Rules for building the application and library
all: build
build:
	$(CC)  $(CFLAGS) $(addprefix src/, $(APP_SRC)) -o $(APP_NAME) $(INCLUDE) $(LIB) -lplplotd -lc -lm -lmpfit

clean:
	@rm -rf $(APP_NAME)
