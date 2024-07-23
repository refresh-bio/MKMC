KMC_DIR_IN = ../kmc

CC = g++
CFLAGS = -Wall -I$(KMC_DIR_IN) -I$(KMC_DIR_IN)/kmcdb/libs -Ilib -O3 -std=c++20
LIBS = $(KMC_DIR_IN)/3rd_party/cloudflare/libz.a
SRCS = $(wildcard *.cpp)
OBJS = $(SRCS:.cpp=.o)

all: mkmc

mkmc: $(OBJS)
		$(CC) $(CLINK) $(OBJS) $(KMC_DIR_IN)/kmc_api/*.o $(KMC_DIR_IN)/kmc_core/*.o $(KMC_DIR_IN)/kmc_dump/nc_utils.o $(LIBS) -o $@

$(OBJS): %.o : %.cpp
		$(CC) $(CFLAGS) -c $< -o $@

clean:
		-rm -f mkmc $(OBJS)
