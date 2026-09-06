CXX = h5cc
LD  = $(CXX)

CXXFLAGS = -O3 -g
LDFLAGS  = -lm $(CXXFLAGS) 

# All .c files in the directory
SRCS := $(wildcard *.c)

# Convert each .c to a .o
OBJS := $(SRCS:.c=.o)

all: eigen_main

# Explicit compile rule for each .c → .o
%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Explicit link rule using all object files
eigen_main: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f eigen_main *.o *~
