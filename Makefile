PLATFORM = $(shell uname)

CFLAGS = -g  -Wall 
LDFLAGS=


ifeq ($(PLATFORM),Darwin)
## Mac OS X
CFLAGS += -m64 -isystem/usr/local/include  -Wno-deprecated 
LDFLAGS+= -m64 -lc -framework AGL -framework OpenGL -framework GLUT -framework Foundation

else
## Linux
CFLAGS += -m64
INCLUDEPATH  = -I/usr/include/GL/ 
LIBPATH = -L/usr/lib64 -L/usr/X11R6/lib
LDFLAGS+=  -lGL -lglut -lrt -lGLU -lX11 -lm  -lXmu -lXext -lXi
endif


CC = g++ -O3 -Wall $(INCLUDEPATH)


PROGS =  main  lidarview
default: $(PROGS)

lidarview: lidarview.o  lidar.o grid.o map.o pixel_buffer.o
	$(CC) -o $@ lidarview.o  lidar.o grid.o map.o pixel_buffer.o  $(LDFLAGS)

main: main.o  lidar.o grid.o map.o pixel_buffer.o
	$(CC) -o $@ main.o lidar.o grid.o map.o pixel_buffer.o  $(LDFLAGS)

lidarview.o: lidarview.cpp lidar.hpp  grid.h map.h pixel_buffer.h 
	$(CC) -c $(INCLUDEPATH) $(CFLAGS)   lidarview.cpp  -o $@

lidar.o: lidar.cpp lidar.hpp grid.h   
	$(CC) -c $(INCLUDEPATH) $(CFLAGS)   lidar.cpp  -o $@

grid.o: grid.c grid.h
	$(CC) -c $(INCLUDEPATH) $(CFLAGS)   grid.c  -o $@

map.o: map.c map.h grid.h pixel_buffer.h stb_image_write.h
	$(CC) $(CFLAGS) -o $@ -c map.c

pixel_buffer.o: pixel_buffer.c pixel_buffer.h stb_image_write.h
	$(CC) $(CFLAGS) -o $@ -c pixel_buffer.c

clean::	
	rm *.o
	rm main 
	rm lidarview

