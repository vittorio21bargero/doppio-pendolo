CXX      = g++
CXXFLAGS ?= -std=gnu++17 -O2 -Wall -fopenmp
INCLUDES  = -Iinclude

SRC = $(wildcard src/*.cpp)
BIN = double_pendulum

$(BIN): $(SRC) $(wildcard include/*.h)
	$(CXX) $(CXXFLAGS) $(INCLUDES) $(SRC) -o $@

run: $(BIN)
	./$(BIN) all

figures:
	python scripts/plot_results.py

clean:
	rm -f $(BIN) $(BIN).exe

.PHONY: run figures clean
