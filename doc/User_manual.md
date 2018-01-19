# CAMUS ver1.0  user's manual

## Environmental setting
  In general, the default stack size of system resources and the default size of the OMP_STACKSIZE environmental variable are as small as 8192 kb.
  Set new stack sizes according to your shell.

  - bash:  
    Write the following two lines in the ~/.bashrc file
    ```
    ulimit -s unlimited  
    export OMP_STACKSIZE=1G  
    ```

  - csh, tcsh:  
    Write the following two lines in the ~/.cshrc or ~/.tcshrc file.  
    ```
    unlimit
    setenv OMP_STACKSIZE 1G
    ```

## How to build CAMUS  
  Execute the following command on the shell.
  ```
  $ make
  ```

## How to run CAMUS
  - OpenMP parallel computing  
    Set the OMP_NUM_THREADS environment variable to specify the number of threads for one process.
  
    - bash: in the ~/.bashrc file
    ```
    export OMP_NUM_THEADS=(# of threads)
    ```

    - csh, tcsh:  in the ~/.cshrc or ~/.tcshrc file
    ```
    setenv OMP_NUM_THEADS (# of threads)
    ```
    
    After the setting, start a calculation.
    ```
    $ bin/camus input > output
    ```
  
##  Input file format
1. input file  
  - input1 namelist  
    - natom          : total number of particle   
    - nstep_total    : total number of simulation steps  
    - dump_ene_step  : energy is dumped at intervals of (dump_ene_step)  
    - dump_step      : structure is dumped at intervals of (dump_step)   
  
  - data_file namelist  
    - pos_file       : position file name  
    - vel_file       : velocity file name  
    - aij_file       : aij file name  
    - particle_file  : particle file name  
    - bond_file      : bond file name  
  
  - phys_param namelist
    - box_length(1) ,box_length(2) ,box_length(3) : length of simulation box in reduced unit
    - dt  : time step in reduced unit
    - kb  : boltzmann constant in reduced unit
    - temp : temperature in reduced unit
    - lambda      : lambda parameter for DPD simulation
  
  - dpd_parameter namelist
    - gamma       : gamma
    - cutoff      : cutoff distance
  
  - calc_flags namelist
    - use_random_initial_coord :  use random structure for initial coordinate
    - temp_control             :  use temperature control ( velocity scaling method )

2. position file  
   file format : particle index , x, y, z
   ```
  # example
  1           7.28882           1.12129           7.37967
   ```

3. velocity file  
   file format : particle index , vx, vy, vz
   ```
  # example
  1          -0.82173          -0.04216           0.41583
   ```

4. aij file  
  line1: comment  
  line2: "particle type1" "particle type2" aij  
   ```
  ## example
  1 1 50.0
  1 2 55.0
  2 2 50.0
   ```

5. particle file  
  particle index, particle type, fix flag ( 0 -> fix, 1 -> movable)  
  line1: comment  
  line2: # particle index, particle type  
   ```
  # example
  1 2 1
  2 2 1
  3 2 1
  4 2 1
   ```

6. bond file  
  bond type (1-> harmonic, 2-> morse), particle index A, particle index B, length, spring const, spring const2 ( morse bond only )  
  line1: comment  
  line2: bond type, particle index A, particle index B, natural length, spring constant, spring constant2 ( morse bond only )  
   ```
  # example
  1   1    2    0.6   160.0   0 <- spring bond
  2   1    3    0.9   12.0    8.0 <- morse bond
   ```

