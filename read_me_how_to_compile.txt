# How to compile the solvers using gcc (should at least be valid for mac and linux)
1) -> open terminal in this folder (RunMCISimulation) or got to this directory via terminal
2) -> run: gcc -shared -fPIC -o output_name.so solver_name.c helpers.c RNG.c
       eg: gcc -shared -fPIC -o linear_solver.so verlet_linear.c helpers.c RNG.c

