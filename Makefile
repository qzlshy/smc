OBJS = main.o readres.o readpdb.o l_math.o readconf.o readtop.o readprm.o build.o caldist.o
main: ${OBJS}
	g++ -o long -lm ${OBJS}

clean:
	rm -f ${OBJS}
