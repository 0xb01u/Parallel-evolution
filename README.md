# Evolution simulation - a parallel approach
This repository, forked from [HylianPablo's original one](https://www.github.com/HylianPablo/Paralela2020), presents a series of different parallel programs optimizing a sequential simulation for cell/bacteria evolution.

These programs were originally develpoed as part of Universidad de Valladolid's course on Parallel Computing (degree on Software Engineeering, academic year 2019-2020), where students would compete in teams of pairs in a contest to obtain the shortest execution time for the simulation (thus obtaining a greater calification). The original programs were developed by **g110 team**.

# Simulation description and optimization techniques
The students were given a sequential C-coded version of the simulation, and their "only" job was to optimize it using 3 different parallel computing paradigms:
 - Shared memory, using **OpenMP**.
 - Distributed memory, using **MPI**.
 - Accelerators, using **CUDA**.

The original description for the sequential program can be found [here](description.md).

# Contest results
 - g110's OpenMP implementation of the simulation was ranked #5 in the OpenMP contest. The contest's test executions can be found [here](/OpenMP/tests.md).
 - g110's MPI implementation of the simulation was ranked #1 in the MPI contest. The contest's test executions can be found [here](/MPI/tests.md).
 - g110's CUDA implementation of the simulation was ranked #1 in the CUDA contest. The contest's test executions can be found [here](/CUDA/tests.md).

# Repository's strucutre
In the root directory are located various folders with different library/platform names. In those, the following files and folders are found:
 - `evolution_<folder_name>_g110.c`: the latest and most optimized implementation of the evolution simulation developed by g110 team, using the libraries/platforms indicated by the folder the file is in.
 - `Makefile`: the makefile to compile the file specified above, using `make`. The output file will be called `evolution`. The makefiles contains a `debug` mode, similar to the makefiles provided to the students; however, the debug option for the evolution simulations will probably not work in any parallel implementation.
 - \*\*\*`g_110_p<some_number>`: The final version of g110's program in the corresponding leaderboard. As the leaderboards had a closing deadline, this file might not be the same as the one described in the first item of this list.
 - \*\*\*`original_src/`: a folder containing the original files provided to the students to develop the corresponding program.
 - \*\*\*`tests.md`: the execution tests for the respective contest.
 - \*\*\*A `pdf` file: g110's report on the corresponding implementation of the simulation, which had to be included to evaluate the assignment.

(Three asterisks denote the file or folder is related only to the original assignment, and thus will not be found in extra implementations of the simulation.)

# Related repositories
[Cerberus](https://www.github.com/0xb01u/Cerberus) is a Discord bot simulating the CUDA's contest's execution queue server, and was used by various teams to program, debug, and test their programs, including g110. The updated list of execution tests used and developed (partly) by g110 to test their parallel programs can be found [in Cerberus' repository](https://github.com/0xb01u/Cerberus/tree/master/tests), as they were integrated within the bot.
