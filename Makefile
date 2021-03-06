# from http://stackoverflow.com/questions/2908057/can-i-compile-all-cpp-files-in-src-to-os-in-obj-then-link-to-binary-in

CPP_FILES := $(wildcard src/*.cpp)
OBJ_FILES := $(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))
LD_FLAGS := -std=c++11 -fopenmp
CC_FLAGS := -std=c++11 -fopenmp
OBJ := obj

tdca: $(OBJ_FILES)
	g++ $(LD_FLAGS) -o $@ $^

obj/%.o: src/%.cpp
	@mkdir -p $(OBJ)
	g++ $(CC_FLAGS) -c -o $@ $<
