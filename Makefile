all:
	gcc -o bloptest -Wall -Ilibpll/include -Llibpll/lib src/bloptest.c -lpll-sse3 -lm -lnlopt
