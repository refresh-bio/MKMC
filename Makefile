all: mkmc

# *** REFRESH makefile utils
include refresh.mk

$(call INIT_SUBMODULES)
$(call INIT_GLOBALS)
$(call CHECK_OS_ARCH, $(PLATFORM))

# *** Project directories
$(call SET_SRC_OBJ_BIN,mkmc,obj,bin)
3RD_PARTY_DIR := ./3rd_party

# *** Project configuration
#$(call CHECK_NASM)
$(call ADD_KMC_LIB, $(3RD_PARTY_DIR)/kmc)
$(call ADD_ZLIB_NG_AS_ZLIB, $(3RD_PARTY_DIR)/zlib-ng-compat)
#$(call PROPOSE_ISAL, $(3RD_PARTY_DIR)/isa-l)
$(call ADD_MIMALLOC, $(3RD_PARTY_DIR)/mimalloc)
#$(call CHOOSE_GZIP_DECOMPRESSION)
$(call ADD_REFRESH_LIB, $(3RD_PARTY_DIR))
$(call ADD_STATS_LIB, $(3RD_PARTY_DIR)/stats)
$(call ADD_ANNOY_LIB, $(3RD_PARTY_DIR)/annoy)
$(call ADD_HNSWLIB_LIB, $(3RD_PARTY_DIR)/hnswlib)
$(call ADD_UMAPPP_LIB, $(3RD_PARTY_DIR)/umappp)
$(call ADD_CPPIRLBA_LIB, $(3RD_PARTY_DIR)/CppIrlba)
$(call ADD_CPPKMEANS_LIB, $(3RD_PARTY_DIR)/CppKmeans)
$(call ADD_AARAND_LIB, $(3RD_PARTY_DIR)/aarand)
$(call ADD_KNNCOLLE_LIB, $(3RD_PARTY_DIR)/knncolle)
$(call ADD_EIGEN_LIB, $(3RD_PARTY_DIR)/eigen)

$(call SET_STATIC, $(STATIC_LINK))
$(call SET_C_CPP_STANDARDS, c11, c++20)
$(call SET_GIT_COMMIT)

$(call SET_FLAGS, $(TYPE))

$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_x86_64, 10, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Linux_aarch64, 11, 20)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_x86_64, 11, 13)
$(call SET_COMPILER_VERSION_ALLOWED, GCC, Darwin_arm64, 11, 13)

ifneq ($(MAKECMDGOALS),clean)
$(call CHECK_COMPILER_VERSION)
endif

# *** Source files and rules
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,MAIN,))
$(eval $(call PREPARE_DEFAULT_COMPILE_RULE,KMC_API,kmc_api))

# *** Targets
mkmc: $(OUT_BIN_DIR)/mkmc
$(OUT_BIN_DIR)/mkmc: zlib-ng libkmc mimalloc_obj \
	$(OBJ_MAIN) $(OBJ_KMC_API)
	-mkdir -p $(OUT_BIN_DIR)	
	$(CXX) -o $@  \
	$(OBJ_MAIN) $(OBJ_KMC_API) \
	$(LIBRARY_FILES) $(LINKER_FLAGS) $(LINKER_DIRS)

# *** Cleaning
.PHONY: clean init
clean: clean-zlib-ng clean-mimalloc_obj
	-rm -r $(OBJ_DIR)
	-rm -r $(OUT_BIN_DIR)

init:
	$(call INIT_SUBMODULES)
