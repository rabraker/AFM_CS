# ROOT=/home/arnold/sandbox/c_sandbox/plplot_sandbox
TARGET = linux
CC  = gcc
CFLAGS= -g -pedantic -Wall -Wextra -Wundef -Winit-self -fPIC

#include makedefs
LIB_SRC=   SOS.c fit_ZPK.c

SRC_DIR = src
BUILD_DIR = build
APP = bode_fit
LIBNAME_ = fit_zpk
CPPFLAGS  = -Icmpfit/ -Iinclude
LDFLAGS= -Lcmpfit/

ifeq ($(TARGET), linux)
	CPPFLAGS += -I/usr/include/plplot/
	LDFLAGS  += -L/usr/lib/x86_64-linux-gnu \
				-L/usr/lib/x86_64-linux-gnu/plplot5.10.0
	LDLIBS = -lplplotd -lc -lm -lmpfit
	LIBNAME = lib$(LIBNAME_).so
	SRC = $(LIB_SRC)
	APP_SRC = src/main.c src/test_data.c src/plotting.c
else
	APP_SRC = src/main.c src/test_data.c
	LD = gcc
	CPPFLAGS += -DBUILDING_DLL -D__MINGW32__
	LDLIBS = -lc -lm -lmpfit
	LIBNAME = lib$(LIBNAME_).dll
endif

OBJ       = $(SRC:%.c=$(BUILD_DIR)/%.o)

LIBPATH = $(shell pwd)


$(info $(OBJ))
# Rules for building the application and library
all: $(LIBNAME) $(APP)


# 	$(CC) $(CPPFLAGS) $(CFLAGS) $(addprefix src/, $(APP_SRC)) -o $(APP_NAME)  $(LIB) -lplplotd -lc -lm -lmpfit

$(APP):$(APP_SRC)
	$(CC) $(CPPFLAGS) $(CFLAGS) -Wl,--rpath $(LIBPATH) $(LDFLAGS) \
		-L./ -o $(APP) $(APP_SRC) -l$(LIBNAME_) $(LDLIBS)

$(LIBNAME):$(OBJ)
	$(CC) -shared $(LDFLAGS) $(OBJ) -o $(LIBNAME) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

clean:
	@rm -f $(APP)
	rm -f $(BUILD_DIR)/*.o
	rm -f $(LIBNAME)
