all: mkmc


dummy := $(shell git submodule update --init --recursive)


MKMC_MAIN_DIR = mkmc
OUT_BIN_DIR = bin

KMC_DIR = 3rd_party/kmc
KMC_ZLIB_DIR = $(KMC_DIR)/3rd_party/cloudflare

KMC_LIB_ZLIB = $(KMC_ZLIB_DIR)/libz.a
LIB_KMC = $(KMC_DIR)/bin/libkmc_core.a

KMC_LIB_NC_UTILS = $(KMC_DIR)/kmc_dump/nc_utils.o


LIBS=-I$(MKMC_MAIN_DIR)/lib \
     -I$(MKMC_MAIN_DIR)/lib/stats/include \
     -I$(MKMC_MAIN_DIR)/lib/annoy/include \
     -I$(MKMC_MAIN_DIR)/lib/hnswlib \
     -I$(MKMC_MAIN_DIR)/lib/umappp/include \
     -I$(MKMC_MAIN_DIR)/lib/CppIrlba/include \
     -I$(MKMC_MAIN_DIR)/lib/CppKmeans/include \
     -I$(MKMC_MAIN_DIR)/lib/aarand/include \
     -I$(MKMC_MAIN_DIR)/lib/knncolle/include \
     -I$(MKMC_MAIN_DIR)/lib/eigen


ifdef MSVC     # Avoid the MingW/Cygwin sections
    UNAME_S := Windows
else                          # If uname not available => 'not'
    UNAME_S := $(shell sh -c 'uname -s 2>/dev/null || echo not')
    UNAME_M := $(shell uname -m)
endif

D_OS =
D_ARCH =

ifeq ($(UNAME_S),Darwin)
	D_OS=MACOS
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	else
		D_ARCH=X64
	endif
else
	D_OS=LINUX
	D_ARCH=X64
	ifeq ($(UNAME_M),arm64)
		D_ARCH=ARM64
	endif
	ifeq ($(UNAME_M),aarch64)
		D_ARCH=ARM64
	endif
endif

CPU_FLAGS =
STATIC_LFLAGS =
PLATFORM_SPECIFIC_FLAGS =

#in some cases we can have different results on ARM
#I guess this is exactly the same as here: https://bugs.mysql.com/bug.php?id=82760
ifeq ($(D_ARCH),ARM64)
	PLATFORM_SPECIFIC_FLAGS = -ffp-contract=off
endif

ifeq ($(D_OS),MACOS)
	CC = g++-11

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8.4-a
	else
		CPU_FLAGS = -m64
	endif
	STATIC_LFLAGS = -static-libgcc -static-libstdc++ -pthread
else
	CC 	= g++

	ifeq ($(D_ARCH),ARM64)
		CPU_FLAGS = -march=armv8-a
		STATIC_LFLAGS = -static-libgcc -static-libstdc++ -lpthread
	else
		CPU_FLAGS = -m64
		STATIC_LFLAGS = -static -Wl,--whole-archive -lpthread -Wl,--no-whole-archive
	endif
endif

CLINK_FABI_VERSION = 
ifeq ($(UNAME_S),Linux)
	CLINK_FABI_VERSION = -fabi-version=6
endif


CFLAGS = -fPIC -Wall -O3 $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++20 -pthread $(LIBS) -I $(KMC_DIR) -fpermissive
CLINK = -lm -lpthread

release: CFLAGS += -DNDEBUG
release: CLINK += $(STATIC_LFLAGS)
release: all

debug: CFLAGS = -fPIC -Wall -O0 -g $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++20 -pthread $(LIBS) -I $(KMC_DIR) -fpermissive
debug: all

CLINK += $(CLINK_FABI_VERSION)


MKMC_SRCS = $(wildcard $(MKMC_MAIN_DIR)/*.cpp)	
MKMC_OBJS = $(MKMC_SRCS:.cpp=.o)

$(MKMC_OBJS): %.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(LIB_KMC):
	cd $(KMC_DIR); $(MAKE) bin/libkmc_core.a

$(KMC_LIB_NC_UTILS): %.o : %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

$(KMC_LIB_ZLIB):
	cd $(KMC_ZLIB_DIR); ./configure; $(MAKE) libz.a

mkmc: $(MKMC_OBJS) $(LIB_KMC) $(KMC_LIB_NC_UTILS) $(KMC_LIB_ZLIB)
	-mkdir -p $(OUT_BIN_DIR)
	$(CC) $(CLINK) $(MKMC_OBJS) $(LIB_KMC) $(KMC_LIB_NC_UTILS) $(KMC_LIB_ZLIB) -o $(OUT_BIN_DIR)/$@

install: all
	install bin/* /usr/local/bin

uninstall:
	-rm -f /usr/local/bin/mkmc

clean:
	-rm -rf $(OUT_BIN_DIR)
	-rm -rf $(KMC_LIB_NC_UTILS)
	-rm -rf $(MKMC_OBJS)
	cd $(KMC_ZLIB_DIR) && $(MAKE) clean
	cd $(KMC_DIR) && $(MAKE) clean
