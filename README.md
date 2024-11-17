# Parallel Karnaugh-Maps

## Team
- **Emily Allendorf (eallendo)**
- **Andrew Wang (herongw)**  

## URL
[https://github.com/KingCashVamp/K-Map-Parallelism](https://github.com/KingCashVamp/K-Map-Parallelism)

## SUMMARY
We are going to implement a Karnaugh map solver to minimize Boolean expressions in digital logic design. We will be performing our algorithm on both GPU and multi-core CPU using language frameworks, CUDA and OpenMP. We will then perform a detailed analysis of both systems’ performance characteristics.

## BACKGROUND
### Importance of Boolean Expression Reduction
Boolean expression reduction is a fundamental process in the design and optimization of digital circuits. Simplifying these expressions minimizes the number of gates required, reducing hardware complexity, power consumption, and overall cost. It is critical for hardware description compilation and FPGA Synthesis where physical resource constraints make optimization vital.

### Methods of Boolean Reduction
#### The Quine-McCluskey Algorithm
The Quine-McCluskey method involves finding prime implicants and determining the minimum set of prime implicants that cover all terms (the minimum that can represent equivalent logic as the input). Its tree-like nature makes it a candidate for parallelism.

#### Karnaugh Maps (K-Maps)
Karnaugh Maps (K-maps) are a visual method for simplifying Boolean functions. They represent the truth table of a function in a grid format, with adjacent cells differing by only one variable. Groups of adjacent 1s (or 0s for the Product of Sums approach instead of the Sum of Products approach) are identified and used to generate simplified expressions. K-maps are intuitive and when applied correctly, provide a minimal Boolean representation. However, K-maps are hard to scale and less common in software. Parallelizing a K-map algorithm, because of its complexity especially for large truth tables, is a non-trivial task worth exploring to see if K-maps could be viably used in boolean reduction.

### Objectives
This project aims to focus on the non-trivial task of parallelizing K-maps, exploring the potential of a traditionally manual approach in a modern computational context. We will develop an algorithmic approach for Boolean expression reduction using K-maps.

#### Pseudo-algorithm:
1. Express the Boolean function: write a starting expression of a truth table to be optimized.
2. Set up the K-map: Arrange the truth table logic into the gridlike structure of a K-map.
3. Identify Adjacency: Locate 1s or 0s, and Xs (don’t care values that are common and can be exploited) that may be grouped together.
4. Identify Optimal Groupings: Find the largest possible groupings that can be made of adjacent entries while following grouping constraints (has to be a group the size of a power of 2).

## THE CHALLENGE
### Dependencies:
- Each thread will be responsible for a region of the K-map so synchronization is required to ensure correctness.
- Groups of 1s on the K-map can overlap which introduces the challenge of communication between the processors when considering all possibilities of an optimal solution (this also creates challenges in splitting up the work among the threads).
- All the threads have to be done computing one solution before moving on to look for the next solution which introduces another synchronization barrier.

### Memory Access Characteristics:
- **Spatial Locality**: The K-map will be stored in a binary matrix and adjacent cells are frequently accessed together when looked at as a group which provides good spatial locality.
- **Temporal Locality**: There will be good temporal locality as well since when looking for the optimal collusion, groupings will be revisited to resolve overlapping or nested groups after each adjustment.

### Communication-to-Computation Ratio:
The ratio depends on the size of the K-map and the number of threads. With larger K-maps, the “don't care” conditions increase communication due to the need to merge group results which can lead to higher communication overhead relative to the computational simplicity of checking individual cells. Communication-to-computation will never be low.

### Divergent Execution:
There will be branching conditions that can lead to divergent execution. When determining whether a cell belongs in a group, if it overlaps with another group, there will be multiple ways of how “don’t care” cells and other 1’s should be grouped together. Each of these ways leads down another path to a new solution.

### Constraints:
- **Irregular Data Access**: The workload may lead to irregular memory access patterns due to how “don’t care” conditions are distributed unevenly which would reduce cache efficiency.
- **Synchronization Overhead**: Every thread is responsible for their own region of their K-map. If multiple threads are used and on larger K-maps, overlapping groups would often introduce themselves. This would lead to a lot of careful synchronization to ensure correctness and reduce conflicts. This would lead to large overhead synchronization costs on larger K-maps.
- **Load Balancing**: It will be difficult to find a solution to how to best evenly distribute the work among the processors. We would need to find a way to optimize assigning processors to K-map groupings by the density and arrangement of the “1s” and “don’t care” cells in the K-maps.
- **Scaling**: The computation cost will increase exponentially with the number of variables in the K-maps which would make it hard to scale efficiently without hitting processing limits (sequential).
- **Memory Bandwidth**: For large K-maps, the memory bandwidth required to handle frequent read/write operations might become a bottleneck due to limited cache storage and memory.

## RESOURCES
Our code base will be from scratch whilst following assignment 4 as a reference for creating our binary grid array. We will likely implement the algorithm for solving a K-map ourselves although there might be algorithms out there that we can start with. We both have learned how to solve K-maps by hand, so we would apply the same process in our algorithm whilst keeping in mind ways in which we could exploit parallelism.

## GOALS AND DELIVERABLES
- **What we plan to achieve**: Implement a K-map solver and parallelize it and demonstrate speedup. Since there are multiple different steps in solving the K-map that may be parallelized, we would like to explore two different approaches (CUDA vs OpenMP) that would target two different parts of the problem.
- **What we hope to achieve**: Compare our program to the widely used Quine-McCluskey algorithm and to a parallel version of this algorithm as well.
- **What we will do if work goes more slowly than expected**: We will limit our K-map solver to only one parallel implementation that shows the most theoretical promise. We will also stick to just one of the two programming frameworks (CUDA or OpenMP) and focus on completing our goals rather than comparing the two.

### Demo Details:
- Our poster will include visual representations of K-maps themselves and their groupings as their visual nature is the reason they are popularly used in education.
- We will include drawings of 2-dimensional K-maps as well as higher-dimensional K-maps that are “stacked” on top of each other.
- Apart from that, we will demo the speedup graphs we generate for our analysis to discuss our findings.

### Analysis Goals:
We aim to explore the potential speedup of a parallel K-map algorithm. In our analysis, we will explore other metrics such as cache misses and workload distribution that may impact our speedup. We would also like to place an emphasis on larger truth tables/higher-dimensional K-maps because this seems to be the main downfall of the K-map approach. If time permits we would also like to compare and contrast the potential of the K-map and Quine-McClusky algorithms.

## PLATFORM CHOICE
We decided to work with CUDA and OpenMP. We will be using the NVIDIA GPU gates machines for CUDA (between 48 and 87).

- **CUDA** allows for massive parallelism which is particularly beneficial for larger K-maps with a lot of individual cells. Since each grid cell needs to be evaluated independently, CUDA is a good fit. The GPU architecture is designed to handle heavy workloads which would allow it to have better performance and scalability.
- **OpenMP** allows us to provide high-level directives to parallelize an algorithm. It is easier to go from our code executing sequentially to executing in parallel. It has built-in support for reduction and critical sections which can help in managing dependencies, such as combining overlapping groups in the K-maps. OpenMP is ideal for workloads with moderate parallelism.

Our goal is to work with them separately and compare early on which one is most optimal for our project’s goals.

## SCHEDULE
**Week 1 & 2**:
- Finish sequential implementation by Friday 11/22.
- OpenMP implementation and milestone report by Wednesday 11/27.
- Begin CUDA implementation over Thanksgiving break.

**Week 3 & 4**:
- Finalize CUDA implementation by Wed 12/4.
- Collect speedup data and start analysis Tues 12/10.
- Make the poster by 12/13.
- Compare to Quine-McClusky if time permits.
- Final touches to report by 12/15.
