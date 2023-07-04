INC=-I/usr/include/
LIBS=-L/usr/lib/x86_64-linux-gnu/

emu.out: emu.c main.c pade.c pk_to_xi.c
	gcc -o emu.out -g main.c emu.c pade.c pk_to_xi.c -I/usr/local/include ${INC} -L/usr/local/lib ${LIBS} -lgsl -lgslcblas -lm 

clean:
	$(RM) emu.out

