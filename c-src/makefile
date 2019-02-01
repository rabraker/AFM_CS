# ROOT=/home/arnold/sandbox/c_sandbox/plplot_sandbox
# TARGET = win32
CC  = gcc
CFLAGS= -g -pedantic -Wall -Wextra -Wundef -Winit-self -fPIC

#include makedefs
LIB_SRC   =   SOS.c fit_ZPK.c mpfit.c

SRC_DIR   = src
BUILD_DIR = build
APP       = bode_fit
LIBNAME_  = fit_zpk
CPPFLAGS  = -Icmpfit/ -Iinclude

LDFLAGS = -L./

ifeq ($(TARGET), linux)
CPPFLAGS += -I/usr/include/plplot/
LDFLAGS  += -L/usr/lib/x86_64-linux-gnu \
			-L/usr/lib/x86_64-linux-gnu/plplot5.10.0
# LDFLAGS  += -Lcmpfit/
LDLIBS = -lplplotd -lc -lm
LIBNAME = lib$(LIBNAME_).so

APP_SRC = src/main.c src/test_data.c src/plotting.c
endif
ifeq ($(TARGET), win32)
APP_SRC   = src/main.c src/test_data.c
LD        = gcc
CPPFLAGS += -DBUILDING_DLL -D__MINGW32__
LDLIBS    = -lm
# -lmpfit
LIBNAME   = lib$(LIBNAME_).dll
# LDFLAGS   += -Lcmpfit-32/
endif
ifeq ($(TARGET), win64)
APP_SRC   = src/main.c src/test_data.c
LD        = gcc
CPPFLAGS +=
LDLIBS    = -lm
LIBNAME   = lib$(LIBNAME_).dll
# LDFLAGS   += -Lcmpfit/
endif
$(info "BUILDING target = $(TARGET)")

OBJ       = $(LIB_SRC:%.c=$(BUILD_DIR)/%.o)

LIBPATH = $(shell pwd)


$(info $(OBJ))
# Rules for building the application and library
all: $(LIBNAME) $(APP)


$(APP):$(APP_SRC)
	@echo "----- Building APP------"
	$(CC) $(CPPFLAGS) $(CFLAGS) -Wl,--rpath $(LIBPATH) $(LDFLAGS) \
		 -o $(APP) $(APP_SRC) -l$(LIBNAME_) $(LDLIBS)

$(LIBNAME):$(OBJ)
	@echo "------ linking lib-----------"
	$(CC) -shared $(LDFLAGS) $(OBJ) -o $(LIBNAME) $(LDLIBS)
	#cp lib$(LIBNAME_).dll lib$(LIBNAME_).so

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c
	@echo "----------- Building lib source"
	$(CC) $(CPPFLAGS) $(CFLAGS) -o $@ -c $<

clean:
	@rm -f $(APP)
	rm -f $(BUILD_DIR)/*.o
	rm -f $(LIBNAME_).dll
	rm -f $(LIBNAME_).so
