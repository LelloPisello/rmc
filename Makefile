CFLAGS = -Wall -O2 -D DEBUG

bin/rmc: src/main.c bin/geo.o bin/manifold.o
	@printf "\nCompiling executable source...\n\t"
	clang -o bin/rmc bin/geo.o bin/manifold.o src/main.c $(CFLAGS) 

bin/geo.o: src/geo.c 
	@printf "\nCompiling source file geo.c -> geo.o...\n\t"
	clang -c -o bin/geo.o src/geo.c $(CFLAGS)

bin/manifold.o: src/manifold.c
	@printf "\nCompiling source file manifold.c -> manifold.c\n\t"
	clang -c -o bin/manifold.o src/manifold.c $(CFLAGS)

.PHONY: run clean

run: bin/rmc 
	@bin/rmc

clean: 
	rm bin/*