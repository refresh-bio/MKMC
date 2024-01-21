KMC_DIR_IN = ../$(KMC_DIR)

CC = g++
CFLAGS = -Wall -I$(KMC_DIR_IN) -Ilibs -O3 -std=c++17
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
