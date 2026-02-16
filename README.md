# pancake
2D numerical solver for thin-film flow.

# Installation
`pancake` is a terminal application written in C++ for Linux. To use it, it requires a C++ compiler, preferably CMake (the guide uses CMake). Requirements are:
- C++ compiler (e.g., GCC, Clang, CMake)
- OpenMP
To compile and run `pancake`, follow these steps:
1. Clone the repository:
   ```bash
   git clone https://github.com/yourusername/pancake.git
   ```
2. Navigate to the project directory:
   ```bash
   cd pancake
   ```
3. Compile the code using CMake:
   ```bash
   mkdir build
   cd build
   cmake ..
   make
   ```

# Usage
To run the application in single process:
```bash
./pancake
```
For parallel runs,
```bash
mpirun -n <number_of_processes> ./pancake
```
The solution is stored in the `results` directory. Open `results/results.pvd` in ParaView for data visualization.
