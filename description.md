(Adapted and translated from the original spanish description.)

# 1. Evolution and natural selection simulation

The students are given a sequential code to simulate in a simplified manner the evolution of a cell/bacteria culture.

In the simulation, a rectangular surface where the culture develops is considered. This is represented as a series of equally spaced culture cells or control spots. A **two dimensional array will be used to represent the quantity of food in each culture cell or spot** at a given time. Another one-dimensional array is used to store the information of a series of **bacterias or cells** that move within the surface, capturing food and dividing when they are sufficiently mature, if they have captured enough food.

The information for a cell includes: its age (simulation steps), its position, the direction it moves in, the quantity of food/energy stored, and some genes, that represent the probability of moving forward, rotating left or rotating right each step of the simulation.

### Spreading food
Each simulation step some food units randomly fall into the surface. The number of food spots and the quantity for each one of them is determined by the program's input arguments. It is possible to define a special region where more food falls than in the rest of the surface.

### Cell movement
Each simulation step the cells decide randomly, guided by their probability genes, whether they advance turning right, left. or in the same direction they had. To advance, they consume 1 unit of energy. If they don't have enough energy, they agonize and eventually die. The ones that survive, share the food located in the culture cell they are in to increase their energy reserves.

### Maturing and dividing
When a cell reaches maturity (a pre-established number of simulation steps), if it has enough energy it can be divided into two equal cells, each one with half the energy. When divided, the two cells' age gets resetted to 1, they advance in different directions, and suffer a genetic mutation. These random mutations modify in a greater or lesser degree the probabilities of moving forward or rotating. After enough simulation steps it is possible to observe the effect of natural selection, where the living cells get specialized to better survive in an artificial environment.

### Examples
Two simulation examples (one without special food spot and another with it) can be found on the [output_example file](output_example.md). The symbols on each control spots have the following meaning:
 - Character `'.'`: Food under 5 units.
 - Character `'*'`: Food between 5 and 10 units.
 - Character `'+'`: Food above 10 units.
 - Number between brackets: `(x)`: Culture cell with x living cells.

# 2. Details in the sequential code

## Program arguments
 - rows, columns: Control spots in the surface (size of the array where the surface's food information is stored).
 - max_iter: Maximum number of simulation iterations to execute.
 - max_food: Maximum food level in any position of the surface. If exceeded, the simulation ends.
 - food_density: Proportion of spots in the surface where each simulation step food falls onto.
 - short_rnd1, 2 and 3: Three _unsigned short_ values used to randomly initializate the random number generators for food spreading and cell features.
 - num_cells: initial number of cells.
 - Optional: arguments for the special region where more food falls: position, size, density and food level.

## Results
The program shows the execution time for the computation part (without the initialization), and a statistic results used to verify the correctness of the simulation.

## Debug mode
If compiled with `-DDEBUG` (or using the correct Makefile mode), a function is activated that, during each simulation step, shows a graphical representation of the situation and some statistic results, as showed previously. It can help detecting the precise moments where errors start appearing, or simply show the evolution of the simulation to understand how it works. At the very end the data for the remaining living cells is shown, so their resulting genes are exposed.

## What to parallelize and optimize
In the code are marked the section of the _main_ function that the students should parallelize and optimize. Any code outside this section used for measurements, initializations, etc., should not be modified. The functions that might be modified or optimized by the students are clearly marked with preceding comments.

The DEBUG functions and code is not used in the contest's simulations, and thus there is no need to worry about them. They are only there to ease the debugging of the modified code.

## Objective
The objective is to parallelize and optimize the code without modifying the applied algorithm nor the simulation's results. It is important to distinguish the code sections that can be directly parallelized, the ones that need modifications, solving race conditions, etc.

## About random number generation and parallelism
The students should be cautious when dealing with random number generation. The generation of a sequence of pseudo-random numbers from a seed or state is inherently sequential. In the program some functions from the C standard libraries are used, which allow to store in a small array the current state of the sequence, used as input to generate the next number. This way multiple generators or independent sequences might be used based on saving their states. Nevertheless, it is critical not to alter the order in which the randomly generated numbers for a sequence are obtained, because the results wouldn't be the same.

## Copyright
(c) 2020 Arturo Gonzalez Escribano
Used in Universidad de Valladolid's course on Parallel Computer, from its Software Engineering degree, academic year 2019-2020.
