python_h_path=$(shell python3 -c 'import os, sysconfig; print(os.path.dirname(sysconfig.get_config_h_filename()))')

Z2: include/Z2.hpp src/Z2.cpp swig/Z2.i
	g++ -Wall -O3 -fPIC -g -c src/Z2.cpp src/special.cpp
	swig -I"$(python_h_path)"  -I./ -Wall -python -c++ swig/Z2.i
	g++ -Wall -O3 -fPIC -I"$(python_h_path)" -c swig/Z2_wrap.cxx 
	g++ -shared -o _Z2.so Z2_wrap.o Z2.o special.o
	mv swig/Z2.py swig/Z2_wrap.cxx .

Z3: include/Z3.hpp src/Z3.cpp swig/Z3.i
	g++ -Wall -O3 -fPIC -g -c src/Z3.cpp src/special.cpp
	swig -I"$(python_h_path)"  -I./ -Wall -python -c++ swig/Z3.i
	g++ -Wall -O3 -fPIC -I"$(python_h_path)" -c swig/Z3_wrap.cxx 
	g++ -shared -o _Z3.so Z3_wrap.o Z3.o special.o
	mv swig/Z3.py swig/Z3_wrap.cxx .

calc-M-eta-to-rho-p-lambda: src/bin/calc-M-eta-to-rho-p-lambda.cpp include/Z2.hpp include/special.hpp src/Z2.cpp src/special.cpp
	g++ -Wall -O3  -fopenmp -o calc-M-eta-to-rho-p-lambda calc-M-eta-to-rho-p-lambda.cpp src/Z2.cpp src/special.cpp -laiw -ldl

rho-p-lambda-f2av: src/bin/rho-p-lambda-f2av.cpp
	g++ -Wall -O3  -fopenmp -o rho-p-lambda-f2av src/bin/rho-p-lambda-f2av.cpp -laiw -ldl

rho-p-lambda-f2av-x: src/bin/rho-p-lambda-f2av-x.cpp
	g++ -Wall -O3  -fopenmp -o rho-p-lambda-f2av-x src/bin/rho-p-lambda-f2av-x.cpp 
