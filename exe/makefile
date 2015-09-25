CC = g++
CCFLAG = -O3
INC = $(PWD)/../include/
RND = $(PWD)/../RNG/

all : worm.exe

worm.exe : worm.cpp $(INC)/* makefile
	$(CC) $(CCFLAG) worm.cpp -o worm.exe -I$(INC) -I$(RND)

clean :
	rm *.exe *.diag *.wld *~ *.info *.backup
