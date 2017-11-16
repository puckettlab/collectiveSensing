EXECUTABLE	:= zonalApp
CCFILES		:= particleSystem.cpp  noise.cpp
# Add cuda source files (compiled with nvcc)
CUFILES		:= particleSystem.cu noise.cu
MODULEDIR := ./modules
LOCAL_PREFIX := /media/data1/aawaz/testing/ultimateCode/ultimate6.1
CUDA_INSTALL_PATH = /usr/local/cuda
PNG_INSTALL_PATH = /media/data1/share/software/build
OBJDIR := ./obj
OBJS +=  $(patsubst %.cpp,$(OBJDIR)/%.cpp.o,$(notdir $(CCFILES)))
OBJS +=  $(patsubst %.c,$(OBJDIR)/%.c.o,$(notdir $(CFILES)))
OBJS +=  $(patsubst %.cu,$(OBJDIR)/%.cu.o,$(notdir $(CUFILES)))
CXX=g++##-4.8 
CC=gcc##-4.8
##CXX := g++
NVCC := $(CUDA_INSTALL_PATH)/bin/nvcc 

OPENCV_LIBS =  `pkg-config --cflags opencv`  -lopencv_core  -lopencv_highgui -lopencv_imgproc

##CXXFLAGS= -O3 -Wall -Wno-deprecated `freetype-config --cflags`
#CXXFLAGS= -Og -Wall -Wno-deprecated `freetype-config --cflags`

CXXFLAGS= -O3 -std=c++11 -pedantic -Wall

NVCCFLAGS := -arch=sm_30 \
-gencode=arch=compute_61,code=sm_61 \

INC= -I. -I$(CUDA_INSTALL_PATH)/include -I$(LOCAL_PREFIX)/include -I$(LOCAL_PREFIX)/kernel -I/usr/local/cuda-8.0/samples/common/inc -I$(PNG_INSTALL_PATH)/include 

LIBS= -L$(LOCAL_PREFIX)/lib/ -L/usr/local/lib $(OPENCV_LIBS) -L$(PNG_INSTALL_PATH)/lib -L$(CUDA_INSTALL_PATH)/lib64  -lcudart -lz -lpng -lpngwriter -lgsl -lgslcblas -lm -lfreetype -lcufft 


all: executable
	@echo $(LOCAL_PREFIX)
executable: ${OBJS} main.cpp 
	@echo $(LOCAL_PREFIX)
	$(CXX) main.cpp -o ${EXECUTABLE} $(CXXFLAGS)  $(INC) ${OBJS}  $(LIBS)

$(OBJDIR)/%.cpp.o: %.cpp
	@echo $(CUDA_INSTALL_PATH)
	$(CXX) -c $(CXXFLAGS) $(INC)   $(LIBS) -o $@ $<

$(OBJDIR)/%.cpp.o : $(MODULEDIR)/%.cpp
	$(CXX) -c $(CXXFLAGS) $(INC) $(LIBS)-o $@ $<

$(OBJDIR)/%.cu.o: %.cu
	$(NVCC) -c $(NVCCFLAGS) $(INC) -o $@ $<

clean:
	rm -f obj/*.o ${EXECUTABLE}

