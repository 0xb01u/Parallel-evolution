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
 - g110's OpenMP implementation of the simulation was ranked #5 in the OpenMP contest.
 - g110's MPI implementation of the simulation was ranked #1 in the MPI contest.
 - g110's CUDA implementation of the simulation was ranked #1 in the CUDA contest.

# Repository's strucutre
WIP.

# Related repositories
[Cerberus](https://www.github.com/0xb01u/Cerberus) is a Discord bot simulating the CUDA's contest's execution queue server, and was used by various teams to program, debug, and test their programs, including g110. The updated list of execution tests used and developed (partly) by g110 to test their parallel programs can be found [in Cerberus' repository](https://github.com/0xb01u/Cerberus/tree/master/tests), as they were integrated within the bot.
