# C compiler
CC = nvcc
CC_FLAGS = -O3 -arch=sm_61

prime: 
	$(CC) $(CC_FLAGS) feyman.cu -o feyman 

clean:
	rm -f feyman
