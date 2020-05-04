# Explanation
This folder contains time measurements made with the programs on the HPC Palma of the University of Münster with a NVIDIA's GeForce RTX 2080 Ti. 

## Musket-Data
For understanding it is essential to mention that run refers to starting the program new and iteration refers to how often as part of one run the route is calculated and the pheromone is calculated.
1. `musket_aggregate_splitkernels.csv` - contains the average of several runs with measurements for distinct parts of the program.
2. `musket_splitkernel.csv` - contains the data for single runs with run times for distinct parts of the program. For all calculations executed iteratively (dependent on the number of iterations) the average is denoted. This program also includes the total runtime. For the total runtime of the paper please refer to the next file. 
3. `musket_total_runtime` - contains the average of several runs of the program. In contrast to the previous program, print statements were excluded since they increase the overall runtime for inessential operations. 
