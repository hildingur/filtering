all: ekf

ekf: dir ekf.o filters.o powell.o f1dim.o linmin.o mnbrak.o brent.o
	g++ bin/ekf.o bin/filters.o bin/powell.o bin/f1dim.o bin/linmin.o bin/mnbrak.o bin/brent.o -o bin/ekf_main

ekf.o: ekf.cpp
	g++ -c ekf.cpp -o bin/ekf.o

filters.o: filters.cpp
	g++ -c filters.cpp -o bin/filters.o

powell.o: recipes/powell.cpp
	g++ -c recipes/powell.cpp -o bin/powell.o

f1dim.o: recipes/f1dim.cpp
	g++ -c recipes/f1dim.cpp -o bin/f1dim.o

linmin.o: recipes/linmin.cpp
	g++ -c recipes/linmin.cpp -o bin/linmin.o

mnbrak.o: recipes/mnbrak.cpp
	g++ -c recipes/mnbrak.cpp -o bin/mnbrak.o

brent.o: recipes/brent.cpp
	g++ -c recipes/brent.cpp -o bin/brent.o

dir:
	mkdir -p bin
clean: 
	rm -v ./bin/*
#ekf:
#	g++ ekf.cpp filters.cpp ./recipes/powell.cpp ./recipes/f1dim.cpp ./recipes/linmin.cpp ./recipes/mnbrak.cpp ./recipes/brent.cpp -o ekf	
