SHELL	 = /bin/tcsh
CMP	 = gfortran -c
LNK	 = gfortran 
OPT	 = -O
LIB_TIME = 


OBJ =  ideal_polytropic_gas.o \
	shock_tube_profiles_i.o \
	exact_Riemann_solver_i.o \
	grid_1d.o


Roe.exe: Roe.o $(OBJ)
	($(LNK) $(OPT) -o Roe.exe Roe.o $(OBJ))

Roe.o: Roe.f90 $(OBJ)  
	$(CMP) $(OPT) $< -o $(*).o


ideal_polytropic_gas.o: ideal_polytropic_gas.f90 
	$(CMP) $(OPT) $< -o $(*).o

shock_tube_profiles_i.o: shock_tube_profiles_i.f90 
	$(CMP) $(OPT) $< -o $(*).o

exact_Riemann_solver_i.o: exact_Riemann_solver_i.f90
	$(CMP) $(OPT) $< -o $(*).o

grid_1d.o: grid_1d.f90
	$(CMP) $(OPT) $< -o $(*).o
