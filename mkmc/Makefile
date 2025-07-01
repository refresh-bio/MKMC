KMC_DIR_IN = ../3rd_party/kmc

CC = g++
CFLAGS = -Wall -I$(KMC_DIR_IN) \
-I$(KMC_DIR_IN)/kmcdb/libs \
-Ilib \
-Ilib/stats/include \
-Ilib/annoy/include \
-Ilib/hnswlib \
-Ilib/umappp/include \
-Ilib/CppIrlba/include \
-Ilib/CppKmeans/include \
-Ilib/aarand/include \
-Ilib/knncolle/include \
-Ilib/eigen \
-O3 -std=c++20

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
