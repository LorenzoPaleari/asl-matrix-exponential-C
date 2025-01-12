# GENERAL CONFIG
OBJ_DIR = obj
SRC_DIR = src
BENCHMARK_DIR = benchmark
BASE_DIR = base
EXECUTABLE = bin
TMP_FILE = tmp.bin
OUTPUT_DIR = data
INIT_SCRIPT = init.sh
TEST_SCRIPT = test_runner.sh

# G++ CONFIG
CXX = g++
CXXFLAGS = -g -D FLOP$(if $(findstring S,$(FLOP)),S) -D ROOFLIN$(if $(findstring E,$(ROOFLIN)),E,_) -D NO_BLA$(if $(findstring S,$(NO_BLA)),S,_) -O3 -ffast-math -march=native -mfma -Wall -Wextra -pedantic -I/opt/intel/oneapi/mkl/2023.0.0/include

# GCC CONFIG
CC = gcc
CFLAGS = -g -D FLOP$(if $(findstring S,$(FLOP)),S) -Wall -O3 -ffast-math -march=native -mfma -mavx2 -mavx -I/opt/intel/oneapi/mkl/2023.0.0/include
LDFLAGS = /opt/intel/oneapi/mkl/2023.0.0/lib/intel64/libmkl_intel_lp64.so.2 /opt/intel/oneapi/mkl/2023.0.0/lib/intel64/libmkl_sequential.so /opt/intel/oneapi/mkl/2023.0.0/lib/intel64/libmkl_core.so -lpthread -lm -ldl

# source and object files, determined by the V argument
V_DIR = $(OBJ_DIR)/$(V)
V_SRC = $(wildcard $(SRC_DIR)/$(V)/*.c)
V_OBJ = $(V_SRC:$(SRC_DIR)/$(V)/%.c=$(V_DIR)/%.o)

# source and object files of the base implementation
BASE_SRC = $(wildcard $(SRC_DIR)/$(BASE_DIR)/*.c)
BASE_OBJ = $(BASE_SRC:$(SRC_DIR)/$(BASE_DIR)/%.c=$(OBJ_DIR)/$(BASE_DIR)/%.o)

# benchmark source and object files
BENCHMARK_SRC = $(wildcard $(SRC_DIR)/$(BENCHMARK_DIR)/*.cpp)
BENCHMARK_OBJ = $(BENCHMARK_SRC:$(SRC_DIR)/benchmark/%.cpp=$(OBJ_DIR)/$(BENCHMARK_DIR)/%.o)

# PYTHON CONFIG
PYTHON_VERSION = 3.9
PYHTON_SCRIPT = solver.py
PYHTON_REQUIREMENTS = requirements.txt
VENVDIR = venv

# RULES
.PHONY: build $(EXECUTABLE)

default: build venv
	./$(EXECUTABLE) $(N) ./$(VENVDIR)/bin/python $(PYHTON_SCRIPT) $(TMP_FILE) $(OUTPUT_DIR)

build: $(EXECUTABLE) init

init:
	. ./$(INIT_SCRIPT)

valgrind: build venv
	valgrind --leak-check=full --show-leak-kinds=all --track-origins=yes ./$(EXECUTABLE) $(N) ./$(VENVDIR)/bin/python $(PYHTON_SCRIPT) $(TMP_FILE) $(OUTPUT_DIR)

roofline: build venv
	advixe-cl -collect=roofline -stacks --enable-cache-simulation -project-dir=MyResults -- ./$(EXECUTABLE) $(N) ./$(VENVDIR)/bin/python $(PYHTON_SCRIPT) $(TMP_FILE) $(OUTPUT_DIR)

$(EXECUTABLE): $(V_OBJ) $(BENCHMARK_OBJ) $(BASE_OBJ)
	echo $(V_OBJ)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

$(V_DIR)/%.o: $(SRC_DIR)/$(V)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/$(BASE_DIR)/%.o: $(SRC_DIR)/$(BASE_DIR)/%.c
	@mkdir -p $(@D)
	$(CC) $(CFLAGS) -c $< -o $@

$(OBJ_DIR)/$(BENCHMARK_DIR)/%.o: $(SRC_DIR)/$(BENCHMARK_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) -c $< -o $@

venv:
	python$(PYTHON_VERSION) -m venv $(VENVDIR)
	./$(VENVDIR)/bin/python -m pip install -r requirements.txt

test:
	./$(TEST_SCRIPT)

clean:
	rm -rf $(EXECUTABLE) $(TMP_FILE)

cleanobj:
	rm -rf $(OBJ_DIR)

cleanall: clean cleanobj
	rm -rf $(VENVDIR) $(OUTPUT_DIR)
