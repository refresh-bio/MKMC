all: mkmc kmc_tools

MKMC_MAIN_DIR = mkmc
ZLIB_DIR = kmc/3rd_party/cloudflare
KMC_DIR = kmc

OUT_BIN_DIR=bin


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


CFLAGS	= -fPIC -Wall -O3 $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(ZLIB_DIR) -I $(KMC_DIR) -fpermissive
CLINK	= -lm -lpthread

release: CLINK = -lm -std=c++17 $(STATIC_LFLAGS)
release: CLINK = -lm -std=c++17 $(STATIC_LFLAGS)

release: CFLAGS	= -fPIC -Wall -O3 -DNDEBUG $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(ZLIB_DIR) -I $(KMC_DIR) -fpermissive
release: all

debug: CFLAGS	= -fPIC -Wall -O0 -g $(PLATFORM_SPECIFIC_FLAGS) $(CPU_FLAGS) -std=c++17 -pthread -I $(ZLIB_DIR) -I $(KMC_DIR) -fpermissive
debug: all

ifeq ($(UNAME_S),Linux)
	CLINK+=-fabi-version=6
endif


LIB_ZLIB=$(ZLIB_DIR)/libz.a
LIB_KMC=$(KMC_DIR)/bin/libkmc_core.a

# default install location (binary placed in the /bin folder)
prefix      = /usr/local

# optional install location
exec_prefix = $(prefix)

$(LIB_ZLIB):
	cd $(ZLIB_DIR); ./configure; make libz.a

$(LIB_KMC):
	cd $(KMC_DIR); $(MAKE) bin/libkmc_core.a

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

kmc_tools: $(OUT_BIN_DIR)/kmc_tools

$(OUT_BIN_DIR)/kmc_tools:
	mkdir -p $(OUT_BIN_DIR)
	(cd $(KMC_DIR); $(MAKE) kmc_tools); cp $(KMC_DIR)/bin/kmc_tools $@

mkmc: $(OUT_BIN_DIR)/mkmc

$(OUT_BIN_DIR)/mkmc:
	mkdir -p $(OUT_BIN_DIR)
	cd $(KMC_DIR); $(MAKE) kmc
	cd $(MKMC_MAIN_DIR) && $(MAKE) KMC_DIR=$(KMC_DIR) CC=$(CC) CLINK=$(CLINK)
	-cp $(MKMC_MAIN_DIR)/mkmc $(OUT_BIN_DIR)


install: all
	install bin/* /usr/local/bin

uninstall:
	-rm /usr/local/bin/mkmc

clean:
	-rm -rf $(OUT_BIN_DIR)
	cd $(MKMC_MAIN_DIR) && $(MAKE) clean
	cd $(KMC_DIR) && $(MAKE) clean
