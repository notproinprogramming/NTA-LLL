CXX      = g++
CXXFLAGS = -std=c++17 -O2 -Wall -Wextra

TARGET  = prog
SRCS    = prog.cpp lll.cpp
OBJS    = $(SRCS:.cpp=.o)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $^ -o $@

%.o: %.cpp lll.hpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: clean
