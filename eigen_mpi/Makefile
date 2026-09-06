CXX = h5pcc
LD  = $(CXX)

CXXFLAGS = -O3
LDFLAGS  = -lm $(CXXFLAGS)

# All .c files in the directory
SRCS := $(wildcard *.c)

# Convert each .c to a .o
OBJS := $(SRCS:.c=.o)

all: eigen_mpi_main

# Explicit compile rule for each .c → .o
%.o: %.c
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Explicit link rule using all object files
eigen_mpi_main: $(OBJS)
	$(LD) -o $@ $(OBJS) $(LDFLAGS)

clean:
	rm -f eigen_mpi_main *.o *~
