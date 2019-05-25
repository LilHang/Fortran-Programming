objects = swms_2d.o input2.o materia2.o orthofem.o \
         output2.o sink2.o solute2.o time2.o watflow2.o

run_swms : $(objects)
	gfortran -g -Wall -o run_swms $(objects)
swms_2d.o : swms_2d.f90
	gfortran -c -Wall swms_2d.f90
input2.o : input2.f90
	gfortran -c  -Wall input2.f90
materia2.o : materia2.f90
	gfortran -c -Wall materia2.f90
orthofem.o : orthofem.f90
	gfortran -c -Wall orthofem.f90
output2.o : output2.f90
	gfortran -c -Wall output2.f90
sink2.o : sink2.f90
	gfortran -c -Wall sink2.f90
solute2.o : solute2.f90
	gfortran -c -Wall solute2.f90
time2.o : time2.f90
	gfortran -c -Wall time2.f90
watflow2.o : watflow2.f90
	gfortran -c -Wall watflow2.f90
.PHONY:clean
clean :
		rm $(objects)

