objects = swms_2d.o input2.o materia2.o orthofem.o \
         output2.o sink2.o solute2.o time2.o watflow2.o

run_swms : $(objects)
	gfortran -g -Wall -o run_swms $(objects)
swms_2d.o : swms_2d.for
	gfortran -c swms_2d.for
input2.o : input2.for
	gfortran -c input2.for
materia2.o : materia2.for
	gfortran -c materia2.for
orthofem.o : orthofem.for
	gfortran -c orthofem.for
output2.o : output2.for
	gfortran -c output2.for
sink2.o : sink2.for
	gfortran -c sink2.for
solute2.o : solute2.for
	gfortran -c solute2.for
time2.o : time2.for
	gfortran -c time2.for
watflow2.o : watflow2.for
	gfortran -c watflow2.for
.PHONY:clean
clean :
		del  run_swms $(objects)

