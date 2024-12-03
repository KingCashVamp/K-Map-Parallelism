# Project Schedule and Details

## Week 1 & 2
- Implemented a sequential version of the Quine-McCluskey algorithm.
- Began work on CUDA and MPI implementations.

## Week 3

### 1st Half
- **By 12/4:**
  - Finish attempt at CUDA implementation.
  - Finish attempt at MPI implementation.
- **By 12/5:**
  - Ensure proper synchronization when finding matches within groups.
  - Assign processors to successive groups to ensure only adjacent groups are compared.
  - Message passing for adjacencies that involve 2 processors.
- **By 12/6:**
  - Ensure proper synchronization after every iteration.
  - Implement ring interconnect message passing to share all group updates.
- **By 12/7:**
  - Finalize CUDA implementation (Andrew).
  - Finalize MPI implementation (Emily).
  - If time permits:
    - Optimizations to prevent idle processors.
    - Address the bottleneck caused by decreasing groups (potentially fewer than the number of processors), leading to workload imbalance. 
    - Case on the number of processors and iterations to improve workload balance.
- **By 12/8:**
  - Collect experimental data for each implementation (Andrew & Emily).

## Week 4
- **By 12/10:**
  - Compare and analyze results of the different implementations (Andrew & Emily).
  - Include this in a draft of our final report.
  - Comparisons to make:
    - Speedup graphs.
    - Cache miss data.
    - Sensitivity analysis on batch size.
- **By 12/11:**
  - Draft a version of the poster (Emily).
- **By 12/13:**
  - Finalize the poster (Andrew).

---

## Work Done So Far
Before the break started, we had made an attempt at writing the sequential version of our K-maps algorithm. However, after receiving feedback on our proposal, we decided to pivot to implementing the Quine-McCluskey algorithm. As a result, we were unable to make as much progress as originally hoped over the holiday break, and we have revised our goals based on the feedback received.

The shift to the Quine-McCluskey algorithm allows us to reduce a Boolean expression with a theoretically infinite number of variables, making our code much more computationally intensive. Additionally, we shifted from developing CUDA and OpenMP implementations to CUDA and MPI implementations to make the design less trivial.

We have completed the sequential version of the algorithm, which includes:
- Parsing.
- Generating the prime implicants.
- Identifying the essential prime implicants.
- Timing to measure computation times.

We have also begun working on our parallel versions using CUDA and MPI. We have split the work so that one of us focuses primarily on CUDA and the other on MPI. While progress has been made in the designs, no implementation has been finalized yet.

---

## Goals and Deliverables
Despite falling behind our original schedule, we believe we will still meet most of our goals and deliverables. Since we are focusing on the Quine-McCluskey algorithm, we have abandoned the comparison between our original K-map algorithm idea and Quine-McCluskey. However, we still plan to compare two different parallel implementations. Below are the updated goals and deliverables:

### Goals
- Implement a parallel version of the Quine-McCluskey algorithm for Boolean expression reduction.
- Develop both a CUDA version and an MPI version and compare them.

### Demo Details
- Create a poster (likely digital) with:
  - Visual representations of the algorithm.
  - Explanation of the problem with K-maps and the process of identifying adjacencies.
- Demo speedup graphs and other generated graphs to showcase analysis and comparison of the two different implementations.

### Analysis Goals
We aim to:
- Explore the potential speedup of a parallel Quine-McCluskey algorithm using MPI and CUDA.
- Assess the scalability of each implementation.
- Analyze workload balance, communication costs, and cache efficiency for each implementation.
- Identify potential bottlenecks.

---

## Issues That Concern Us the Most
- **Resources:** We are not worried about resource availability, as we will use CMU's resources available throughout the course.
- **Time:** Our most pressing concern is managing time effectively to complete the work. We have considered our other commitments and assignments when creating this schedule. With a clearer project plan post-proposal, we are optimistic about accomplishing our outlined goals.
