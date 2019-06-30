SRCS = serial.c for.c reduction.c critical.c 201614790222.c

all:$(SRCS)
	gcc serial.c -o serial
	gcc for.c -o for -fopenmp
	gcc reduction.c -o reduction -fopenmp
	gcc critical.c -o critical -fopenmp
	gcc 201614790222.c -o 201614790222 -fopenmp

serial:
	gcc serial.c -o serial

for:
	gcc for.c -o for -fopenmp

	gcc reduction.c -o reduction -fopenmp

critical:
	gcc critical.c -o critical -fopenmp

201614790222:
	gcc 201614790222.c -o 201614790222 -fopenmp

clean:
	rm -rf serial \
	rm -rf for \
	rm -rf reduction \
	rm -rf critical \
	rm -rf 201614790222
