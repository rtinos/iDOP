iDOP_knap : global.o file_man.o fitness.o idop_knap.o immigrants.o selection.o statistics.o transformation.o util_functions.o
	g++ -Wall global.o file_man.o fitness.o idop_knap.o immigrants.o selection.o statistics.o transformation.o util_functions.o -o iDOP_knap

global.o : global.cpp	
	g++ -Wall -o global.o -c global.cpp

file_man.o : file_man.cpp	
	g++ -Wall -o file_man.o -c file_man.cpp

fitness.o : fitness.cpp	
	g++ -Wall -o fitness.o -c fitness.cpp

idop_knap.o : idop_knap.cpp	
	g++ -Wall -o idop_knap.o -c idop_knap.cpp

immigrants.o : immigrants.cpp	
	g++ -Wall -o immigrants.o -c immigrants.cpp

selection.o : selection.cpp	
	g++ -Wall -o selection.o -c selection.cpp

statistics.o : statistics.cpp	
	g++ -Wall -o statistics.o -c statistics.cpp

transformation.o : transformation.cpp	
	g++ -Wall -o transformation.o -c transformation.cpp

util_functions.o : util_functions.cpp	
	g++ -Wall -o util_functions.o -c util_functions.cpp

