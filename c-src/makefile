# Detect what we are compiling on.

SHELL = bash
UNAME_S = $(shell uname -s)
$(info "uname = $(UNAME_S)")

ifeq ($(shell [[ $(UNAME_S) =~ ^MINGW64+ ]] && echo yes || echo no), yes)
	ARCH=win64
endif

ifeq ($(shell [[ $(UNAME_S) =~ ^MINGW32+ ]] && echo yes || echo no), yes)
	ARCH=win32
endif

arch = $(shell [[ $(UNAME_S) =~ ^Linux+ ]] && echo yes || echo no)
$(info "arch = $(arch)" )


ifeq ($(shell [[ $(UNAME_S) =~ ^Linux+ ]] && echo yes || echo no), yes)
	ARCH=Linux
endif

$(info "target = $(ARCH)")

CC  	  = gcc
LD        = gcc

CFLAGS    = -g -pedantic -Wall -Wextra -Wundef -Winit-self -fPIC

#include makedefs
LIB_SRC   = SOS.c fit_ZPK.c mpfit.c
SRC_DIR   = src
APP_SRC   = main.c test_data.c

APP       = bode_fit
LIBNAME_  = fit_zpk
CPPFLAGS  = -Icmpfit/ -Iinclude

LDFLAGS   = -L./

LDLIBS    = -lm


ifeq ($(ARCH), Linux)
	CPPFLAGS += -I/usr/include/plplot/
	LDFLAGS  += -L/usr/lib/x86_64-linux-gnu \
				-L/usr/lib/x86_64-linux-gnu/plplot5.10.0
	LDLIBS    += -lplplotd -lc
	LIBNAME   = lib$(LIBNAME_).so

	APP_SRC  +=  plotting.c
	BUILD_DIR = build
endif

ifeq ($(ARCH), win32)
	BUILD_DIR = build-w32
	CPPFLAGS += -DBUILDING_DLL -D__MINGW32__
	LIBNAME   = lib32$(LIBNAME_).dll
endif

ifeq ($(ARCH), win64)
	BUILD     = build-w64
	LIBNAME   = lib64$(LIBNAME_).dll
endif

$(info "BUILDING target = $(ARCH)")

OBJ       = $(LIB_SRC:%.c=$(BUILD_DIR)/%.o)
APP_OBJ   = $(APP_SRC:%.c=$(BUILD_DIR)/%.o)

LIBPATH   = $(shell pwd)



# Rules for building the application and library
all: $(LIBNAME) $(APP)

$(BUILD_DIR):
	mkdir $(BUILD_DIR)

$(APP):$(APP_OBJ)
	@echo "----- Building APP------"
	$(CC) $(CPPFLAGS) $(CFLAGS) -Wl,--rpath $(LIBPATH) $(LDFLAGS) \
		 -o $(APP) $(APP_OBJ) -l$(LIBNAME_) $(LDLIBS)

$(LIBNAME):$(OBJ)
	@echo "------ linking lib-----------"
	$(CC) -shared $(LDFLAGS) $(OBJ) -o $(LIBNAME) $(LDLIBS)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c |$(BUILD_DIR)
	@echo "----------- Building lib source"
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<


clean:
	rm -f $(APP)
	rm -f $(BUILD_DIR)/*.o
	rm -f $(LIBNAME_).dll
	rm -f $(LIBNAME_).so
