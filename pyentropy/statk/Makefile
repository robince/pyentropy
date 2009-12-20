# Development makefile only for mac
ARCH=-arch i386
CC=/usr/bin/gcc-4.0
CXX=/usr/bin/g++-4.0

all: wrap

wrap: entropy_nsb_c.o wrap.o entropy_c.o gen_c.o sort_c.o hist_c.o
	${CC} ${ARCH} -isysroot /Developer/SDKs/MacOSX10.4u.sdk -g -bundle -undefined dynamic_lookup -L/opt/local/lib -lgsl sort_c.o gen_c.o entropy_c.o entropy_nsb_c.o hist_c.o wrap.o -o wrap.so

wrap.o: wrap.c
	${CC} ${ARCH} -isysroot /Developer/SDKs/MacOSX10.4u.sdk -fno-strict-aliasing -Wno-long-double -no-cpp-precomp -mno-fused-madd -fno-common -dynamic -DNDEBUG -g -O3 -I/Library/Frameworks/Python.framework/Versions/2.5/lib/python2.5/site-packages/numpy/core/include -I/Library/Frameworks/Python.framework/Versions/2.5/include/python2.5 -c wrap.c -o wrap.o

wrap.c:
	cython wrap.pyx

entropy_nsb_c.o:
	${CXX} ${ARCH} -I/opt/local/include -c entropy_nsb_c.cpp

entropy_c.o:
	${CC} ${ARCH} -c entropy_c.c

gen_c.o:
	${CC} ${ARCH} -c gen_c.c

sort_c.o:
	${CC} ${ARCH} -c sort_c.c

hist_c.o:
	${CC} ${ARCH} -c hist_c.c

clean:
	rm -rf *.o wrap.so

cleanpyx:
	rm wrap.c wrap.o wrap.so


