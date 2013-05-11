all: ekf

ekf: ekf.o filters.o powell.o f1dim.o linmin.o mnbrak.o brent.o
	g++ bin/ekf.o bin/filters.o bin/powell.o bin/f1dim.o fin/linmin.o bin/mnbrak.o bin/brent.o -o bin/ekf

bin/ekf.o: ekf.cpp
	g++ -c ekf.cpp -o bin/ekf.o


bin/filters.o: filters.cpp
	g++ -c filters.cpp -o bin/filters.o

bin/powell.o: recipies/powell.cpp
	g++ -c recipies/powell.cpp -o bin/powell.o






#ekf:
#	g++ ekf.cpp filters.cpp ./recipes/powell.cpp ./recipes/f1dim.cpp ./recipes/linmin.cpp ./recipes/mnbrak.cpp ./recipes/brent.cpp -o ekf	
