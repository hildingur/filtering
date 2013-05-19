#alias the g++ command so that we don't have to repecify
compile=g++ -c -o3

#all: ekf ukf heston_ekf_vary_data_qty_bin

#ekf: dir heston_ekf.o filters.o powell.o f1dim.o linmin.o mnbrak.o brent.o filter_utils.o
#	g++ bin/heston_ekf.o bin/filters.o bin/powell.o bin/f1dim.o bin/linmin.o bin/mnbrak.o bin/brent.o bin/choldc.o bin/filter_utils.o -o bin/heston_ekf

#heston_ekf.o: heston_ekf.cpp
#	${compile} heston_ekf.cpp -o bin/heston_ekf.o

#ukf: dir heston_ukf.o filters.o powell.o f1dim.o linmin.o mnbrak.o brent.o filter_utils.o
#	g++ bin/heston_ukf.o bin/filters.o bin/powell.o bin/f1dim.o bin/linmin.o bin/mnbrak.o bin/brent.o bin/choldc.o bin/filter_utils.o -o bin/heston_ukf

#heston_ekf_vary_data_qty_bin: dir heston_ekf_vary_data_qty.o filters.o powell.o f1dim.o linmin.o mnbrak.o brent.o filter_utils.o
#	    g++ bin/heston_ekf_vary_data_qty.o bin/filters.o bin/powell.o bin/f1dim.o bin/linmin.o bin/mnbrak.o bin/brent.o bin/choldc.o bin/filter_utils.o -o bin/heston_ekf_vary_data_qty

#heston_ukf.o: heston_ukf.cpp
#	${compile} heston_ukf.cpp -o bin/heston_ukf.o

#heston_ekf_vary_data_qty.o: heston_ekf_vary_data_qty.cpp
#	${compile} heston_ekf_vary_data_qty.cpp -o bin/heston_ekf_vary_data_qty.o


all: final_ekf

final_ekf: dir filters.o powell.o filter_utils.o final_ekf.o
	g++ bin/filters.o bin/powell.o bin/f1dim.o bin/linmin.o bin/mnbrak.o bin/brent.o bin/choldc.o bin/filter_utils.o bin/final_ekf.o -o bin/final_ekf

final_ekf.o: final_ekf.cpp
	${compile} final_ekf.cpp -o bin/final_ekf.o

filters.o: filters.cpp
	${compile} filters.cpp -o bin/filters.o

powell.o: f1dim.o linmin.o mnbrak.o brent.o choldc.o recipes/powell.cpp
	${compile} recipes/powell.cpp -o bin/powell.o

f1dim.o: recipes/f1dim.cpp
	${compile} recipes/f1dim.cpp -o bin/f1dim.o

linmin.o: recipes/linmin.cpp
	${compile} recipes/linmin.cpp -o bin/linmin.o

mnbrak.o: recipes/mnbrak.cpp
	${compile} recipes/mnbrak.cpp -o bin/mnbrak.o

brent.o: recipes/brent.cpp
	${compile} recipes/brent.cpp -o bin/brent.o

choldc.o: recipes/choldc.cpp
	${compile} recipes/choldc.cpp -o bin/choldc.o

filter_utils.o: filter_utils.h filter_utils.cpp choldc.o
	${compile} filter_utils.cpp -o bin/filter_utils.o



dir:
	mkdir -p bin
clean: 
	rm -v ./bin/*
#ekf:
#	g++ ekf.cpp filters.cpp ./recipes/powell.cpp ./recipes/f1dim.cpp ./recipes/linmin.cpp ./recipes/mnbrak.cpp ./recipes/brent.cpp -o ekf	
