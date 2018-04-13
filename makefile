SOURCE = ./source
INCLUDE = ./include
LIB = ./lib
DOC = ./doc
RUN = ./run
MAIN = ./main



pgmk:
	g++ -c $(SOURCE)/alignment.cpp $(SOURCE)/graph.cpp $(SOURCE)/rg.cpp $(SOURCE)/edge.cpp -I$(INCLUDE) -std=c++11
	rm -f $(LIB)/libgraph.a
	ar r $(LIB)/libgraph.a alignment.o graph.o rg.o edge.o
	rm *.o
	g++ $(MAIN)/main53.cpp -L$(LIB) -I$(INCLUDE) -lgraph -o $(RUN)/pgmk -std=c++11 -O3



mproper:
	g++ -c $(SOURCE)/alignment.cpp $(SOURCE)/graph.cpp $(SOURCE)/edge.cpp -I$(INCLUDE) -std=c++11
	rm -f $(LIB)/libgraph.a
	ar r $(LIB)/libgraph.a alignment.o graph.o edge.o
	rm *.o
	g++ $(MAIN)/main.cpp -L$(LIB) -I$(INCLUDE) -lgraph -o $(RUN)/mproper -std=c++11 -O3