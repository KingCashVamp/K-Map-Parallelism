# Parallel Karnaugh Maps

### **Team**  
- Emily Allendorf (eallendo)  
- Andrew Wang (herongw)  

### **Project Repository**  
[GitHub Repository](https://github.com/KingCashVamp/K-Map-Parallelism)

---

## **Summary**  
We are implementing a Karnaugh map solver to minimize Boolean expressions in digital logic design. This algorithm will be executed on both GPU and multi-core CPU platforms using CUDA and OpenMP, respectively. We aim to analyze and compare the performance characteristics of these two systems.

---

## **Background**  

### Importance of Boolean Expression Reduction  
Boolean expression reduction is fundamental in the design and optimization of digital circuits. Simplified expressions minimize the number of gates, reducing hardware complexity, power consumption, and overall cost. This process is especially critical in hardware description compilation and FPGA synthesis where physical resource constraints make optimization essential.

### Methods of Boolean Reduction  
#### **The Quine-McCluskey Algorithm**  
This method finds prime implicants and determines the minimum set of implicants that cover all terms. Its tree-like structure makes it suitable for parallelization.  

#### **Karnaugh Maps (K-Maps)**  
K-maps are a visual tool for simplifying Boolean functions. They represent truth tables in a grid format, where adjacent cells differ by only one variable. Groups of adjacent `1s` (or `0s` for POS) are identified to generate simplified expressions. While intuitive, K-maps are difficult to scale and are less common in software. This project explores the non-trivial task of parallelizing K-map algorithms to assess their viability in modern computation.

---

## **Objectives**  
We aim to develop an algorithmic approach for Boolean expression reduction using K-maps.  

### **Pseudo-Algorithm**  
1. **Express the Boolean Function**: Write a starting expression of a truth table to be optimized.  
2. **Set up the K-map**: Arrange the truth table logic into the grid structure of a K-map.  
3. **Identify Adjacency**: Locate adjacent `1s`, `0s`, and `X` (don't care) values.  
4. **Identify Optimal Groupings**: Find the largest adjacent groupings that follow constraints (group sizes must be powers of 2).  

---

## **The Challenge**  

### Dependencies  
- Synchronization is required to ensure correctness as threads handle overlapping regions of the K-map.  
- Synchronization barriers are needed for all threads to complete one solution before starting the next.  

### Memory Access Characteristics  
- **Spatial Locality**: Adjacent cells in the K-map grid are accessed frequently.  
- **Temporal Locality**: Groupings are revisited for resolving overlaps and adjustments.  
- **Communication-to-Computation Ratio**: Increases with larger K-maps due to merging results and handling "don't care" conditions.  

### Constraints  
- **Irregular Data Access**: Uneven distribution of "don't care" conditions may reduce cache efficiency.  
- **Synchronization Overhead**: Frequent communication among threads introduces latency.  
- **Load Balancing**: Assigning processors efficiently to regions of the K-map is challenging.  
- **Scaling**: Computational cost grows exponentially with the number of variables.  
- **Memory Bandwidth**: Frequent read/write operations for large K-maps may bottleneck performance.  

---

## **Resources**  
- Our implementation will be developed from scratch, using prior coursework (Assignment 4) as a reference for creating the binary grid array.  
- Algorithms for solving K-maps will follow the manual process we’ve learned, adapted for parallelism.  

---

## **Goals and Deliverables**  

### Planned Deliverables  
- A functional K-map solver.  
- Parallel implementations using OpenMP and CUDA.  
- Demonstration of speedup.  

### Aspirational Goals  
- Compare performance with the Quine-McCluskey algorithm and its parallel variant.  

### Contingency Plans  
- If progress slows, focus on a single parallel implementation (OpenMP or CUDA).  

---

## **Demo Details**  
Our demo will include:  
- Visual representations of 2D and higher-dimensional K-maps.  
- Speedup graphs and analysis of our findings.  

---

## **Platform Choice**  
- **CUDA**: Optimized for high parallelism in evaluating grid cells independently, benefiting large K-maps.  
- **OpenMP**: Ideal for moderate parallelism with built-in support for reduction and critical sections.

---

## **Schedule**  

### Week 1 & 2 (11/10 – 11/27):  
- Finish sequential implementation by **11/22**.  
- Implement OpenMP version and complete the milestone report by **11/27**.  
- Begin CUDA implementation over Thanksgiving break.  

### Week 3 & 4 (11/28 – 12/15):  
- Finalize CUDA implementation by **12/4**.  
- Collect speedup data and start analysis by **12/10**.  
- Create the poster by **12/13**.  
- Compare with Quine-McCluskey algorithm if time permits.  
- Finalize report by **12/15**.  

