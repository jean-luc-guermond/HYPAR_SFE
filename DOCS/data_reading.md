Data reading and writing
========================


A built-in user data framework is provided with HYPAR_SFE. To build one's own arguments:

### Initializing the program

Besides the initialization of mpi and petsc, the program must initialize the data reading framework so that the first lines it reads should roughly be:

    USE petsc
    USE read_inputs_module, ONLY:  clean_data_once

    !===Start PETSC and MPI (mandatory)
    CALL PetscInitialize(PETSC_NULL_CHARACTER, ierr)
    communicator = PETSC_COMM_WORLD
    CALL MPI_Comm_rank(communicator, rank, ierr)

    !===Clean data once
    CALL clean_data_once

### Create an argument type

Define one's own my_argument_type. It contains all chains of characters that the user wishes to see inside the data file. 
**Warning**: an argument must absolutely start with "=== " (the last space character is mandatory here) and finish with " ===" (the first space character is mandatory)

    USE read_inputs_module, ONLY: rec_length

    TYPE my_argument_type
      CHARACTER(LEN=rec_length) :: arg1 = "=== some argument ===" #is legal 
      CHARACTER(LEN=rec_length) :: arg2 = "===another argument ===" #is illegal 
      CHARACTER(LEN=rec_length) :: arg3 = "=== another argument===" #is illegal 
      CHARACTER(LEN=rec_length) :: arg4 = "===another argument===" #is illegal 
    END TYPE my_argument_type

<div><u>Remark</u></div>rec_length = 200 is the maximum length of an argument. If required rec_length can be increased manually within HYPAR_SFE/COMMON_SOURCES/TOOLS/read_inputs_module.F90

### Create a user_parameter type

Create a second type which will receive the information read by the user. At the declaration point, it is advised to initialize already the variables. It attributes a default value to the variable which will be used if the corresponding argument line isn't found within the data file

**WARNING**: to read a chain of characters from the data file, it must be declared as

        CHARACTER(LEN=rec_length) :: some_variable

More generally, the types that can be read can be found in "INTERFACE read_data" (HYPAR_SFE/COMMON_SOURCES/TOOLS/read_inputs_module.F90). An example of a parameter type is:

    USE read_inputs_module, ONLY: rec_length

    TYPE my_argument_type
      CHARACTER(LEN=rec_length) :: val1 = "=== 1st argument ===" 
      CHARACTER(LEN=rec_length) :: val2 = "=== 2nd argument ===" 
      CHARACTER(LEN=rec_length) :: val3 = "=== 3rd argument ===" 
      CHARACTER(LEN=rec_length) :: val4 = "=== 4th argument ===" 
    END TYPE my_argument_type

    TYPE my_parameter_type
      CHARACTER(LEN=rec_length)               :: val1 = "method_1" 
      REAL(KIND=8)                            :: val2 = 1.d0
      INTEGER                                 :: val3 = 3 
      REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: val4 !should be allocated with correct size before reading
    CONTAINS
      PROCEDURE :: read => read_my_data
    END TYPE my_parameter_type

### Write the subroutine which will do the reading

It needs to have this precise structure:

    SUBROUTINE read_my_data(this)
        USE read_inputs_module, ONLY: read_data_init_list, read_data, finalize_rewrite_data
        IMPLICIT NONE
        CLASS(my_parameter_type), INTENT(INOUT) :: this
        TYPE(my_argument_type)                  :: argument_data
    
        !=== Reading all data file
        !=== section_name is optional, and only serves for data writing
        CALL read_data_init_list(section_name)

        !=== reading some arguments
        CALL read_data(argument_data%val1, this%val1)
        CALL read_data(argument_data%val2, this%val2, opt_add=this%val1==1)
        CALL read_data(argument_data%val3, this%val3, opt_name='euler')
        ALLOCATE(this%var4(this%val3))
        CALL read_data(argument_data%val4, this%val4, opt_name='euler')

        !=== Closes/Rewrites the data file with what was read
        CALL finalize_rewrite_data

    END SUBROUTINE read_my_data

1. Start with **CALL read_data_init_list** and end with **CALL finalize_rewrite_data**. 

2. In *read_data_init_list*, *section_name* is optional and will define a section in your data file.

3. Repeat the calls **CALL read_data(argument_data%some_parameter, this%some_parameter)**.

4. By default, if the line *argument_data%some_parameter* is absent from the data file being currently read, it is added automatically. This default behavior is disabled by setting the optional *opt_add=some_logical_condition*. The same purpose can be achieved with *if* blocks, but this solution is strongly discouraged since it may completely rearrange the data file in an undesirable way.

5. The user might want to add the same line *argument_data%some_parameter* several times but with some tiny variants. The optional *opt_name* exists for that purpose: the line *CALL read_data(=== some_argument ===, some_val, opt_name="with some additional string")* will actually read the line *=== some_argument with some additional string ===*. An example of usage is for Dirichlet boundary conditions, which may be specified for different fields.

6. The data files that are read/wrote can be through srun -n1 ./a.exe data_in=some_file data_out=another_file data_save=a_third_file

> data_in: the file read by the program. data_in=data by default

> data_out: the modified file wrote by the program. data_out=data by default, meaning data is overwritten by default.

> data_save: a copy of data_in. data_save=previous_data by default. It is relevant if data_in=data_out and in case the code runs into a bug while reading the data file.

Non-exhaustive list of warnings
===============================

1. If the program somehow fails during the data reading, the new data file may be badly rearranged. A previous version of the data file can therefore be found within "previous_data" **FIXME** which is overwritten at every run.

2. Before doing any data reading, the program must absolutely call **CALL clean_data_once**. This call should always be done right *after* initializing mpi and petsc.

That's it! The simplest example can be found in HYPAR_SFE/COMMON_SOURCES/SOLVERS/solver_data_module.F90. Most applications have their own *start_setup.F90* module which itself does some data reading and can be completed with one's own requirements.