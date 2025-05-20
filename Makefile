CXX = g++
#CXXFLAGS = -v -std=c++11 -O3 -march=native -ftree-vectorize
CXXFLAGS = -std=c++11 -O3 -march=native -ftree-vectorize
DEBUGFLAGS = -std=c++11 -g -Og

OBJS = run_mgs.o mulguisin_kd.o
TARGET = run_mgs

INCLUDE = -Iqhull/include -I/home/young/hdf5/include
LDFLAGS = -lqhullcpp -lqhullstatic_r -lhdf5_cpp -lhdf5
LIBS = -Lqhull/lib -L/home/young/hdf5/lib

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $^ $(LIBS) $(LDFLAGS)

%.o: src/%.cc
	$(CXX) $(CXXFLAGS) $(INCLUDE) -c $< -o $@

clean:
	rm -f *.o $(TARGET)
