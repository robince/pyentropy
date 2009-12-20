
ARCH=-arch i386
CC=/usr/bin/gcc-4.0
CXX=/usr/bin/g++-4.0

all: nsb

nsb: entropy_nsb.o nsb.o
	${CC} ${ARCH} -isysroot /Developer/SDKs/MacOSX10.4u.sdk -g -bundle -undefined dynamic_lookup -L/opt/local/lib -lgsl sort_c.o gen_c.o entropy_c.o nsb_entropy.o nsb.o -o nsb.so

nsb.o: nsb.c
	${CC} ${ARCH} -isysroot /Developer/SDKs/MacOSX10.4u.sdk -fno-strict-aliasing -Wno-long-double -no-cpp-precomp -mno-fused-madd -fno-common -dynamic -DNDEBUG -g -O3 -I/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy/core/include -I/Library/Frameworks/Python.framework/Versions/2.5/include/python2.5 -c nsb.c -o nsb.o

nsb.c:
	cython nsb.pyx

entropy_nsb.o: nsb_entropy.o entropy_c.o gen_c.o sort_c.o
	#${CXX} ${ARCH} sort_c.o gen_c.o entropy_c.o nsb_entropy.o -shared -L/opt/local/lib -lgsl -o entropy_nsb.o

nsb_entropy.o:
	${CXX} ${ARCH} -I/opt/local/include -c nsb_entropy.cc

entropy_c.o:
	${CC} ${ARCH} -c entropy_c.c

gen_c.o:
	${CC} ${ARCH} -c gen_c.c

sort_c.o:
	${CC} ${ARCH} -c sort_c.c

clean:
	rm -rf *.o nsb.c nsb.so


