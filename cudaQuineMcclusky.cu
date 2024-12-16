#include <string>
#include <algorithm>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <vector>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <cstring>
#include <unistd.h>
#include <omp.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <driver_functions.h>

/**
 * CUDA K-Map Reduction via QuineMcCluskey algorithm
 * Andrew Wang(herongw)
 */

const int MAX_TERMS = 1024; // Maximum number of terms allowed
const int MAX_VARS = 11;   // Maximum number of variables allowed
#define BLOCKSIZE 1024

struct term {
    char stringTerm[MAX_VARS + 1]; // +1 for null terminator
    int decimals[MAX_TERMS];
    int decimalCount;
    int tickmark;
    int groupNum;
};

struct group {
    term terms[MAX_TERMS];
    int size; //the number of terms in the group
};

struct Lock {
    int* mutex;

    Lock() {
        int state = 0;
        cudaMalloc((void**)&mutex, sizeof(int));
        cudaMemset(mutex, 0, sizeof(int)); // Initialize mutex to 0
        cudaMemcpy(mutex, &state, sizeof(int), cudaMemcpyHostToDevice);
    }

    ~Lock() {
        cudaFree(mutex);
    }

    __device__ void lock() {
        while (atomicCAS(mutex, 0, 1) != 0) {
            // Busy-wait until lock is acquired
        }
    }

    __device__ void unlock() {
        atomicExch(mutex, 0);
    }
};

//very basic command line parser to process the given information
// example command line for now: ./Kmap n 4 t 1 2 3 15 16 ...
int parseInput(int argc, char* argv[], std::vector<int>& givens){

    // printf("0: %s, 1: %s , 2: %s\n", argv[0], argv[1], argv[2]);
    int i = 1;
    int num_vars;
    if (strcmp(argv[i],"n") == 0) i++;
    else {
        printf("error 0: parsing the command line\n");
        return -1;
    }
    
    try {
        num_vars = std::stoi(argv[i]);
        i++;
        } catch (const std::invalid_argument& error) {
            printf("error 1: parsing the command line\n");
            return -1;
        } catch (const std::out_of_range& error){
            printf("error 2: parsing the command line\n");
            return -1;
        }
    
    if (strcmp(argv[i],"t") == 0) i++;
    else {
        printf("error 3: parsing the command line\n");
        return -1;
    }
    // printf("GOTHERE");
    for (int j = i; j < argc; j++) {
       try {
            int term = std::stoi(argv[j]);
            givens.push_back(term);
        } catch (const std::invalid_argument& error) {
            printf("error 4: parsing the command line\n");
            return -1;
        } catch (const std::out_of_range& error){
            printf("error 5: parsing the command line\n");
            return -1;
        }
    }

    return num_vars;
}

__device__  term setMatch(term a, term b, int num_vars, int* difference) {
    //int index = blockIdx.x * blockDim.x + threadIdx.x;
    term res;

    int ones = 0;
    res.decimalCount = 0;
    
    for (int i = 0; i < num_vars; i++){
        if (a.stringTerm[i] != b.stringTerm[i]) {
            (*difference)++;
            res.stringTerm[i] = '_';
        } else {
            const char bin = a.stringTerm[i];
            res.stringTerm[i] = bin;
            if (bin == '1') ones += 1;
        }
        res.decimalCount++;
    }
    res.stringTerm[num_vars] = '\0'; // Null-terminate the string
    
    //populate the result decimals with all of the minterms involved
    for (int i = 0; i < a.decimalCount; i++) {
        res.decimals[res.decimalCount] = a.decimals[i];
        res.decimalCount++;
    }
    for (int i = 0; i < b.decimalCount; i++) {
        res.decimals[res.decimalCount] = b.decimals[i];
        res.decimalCount++;
    }
  

    res.tickmark = 0;
    res.groupNum = ones;

    // if (index == 0) {
    //     printf("A string =");
    //     for (int i = 0; i < num_vars; i++) {
    //         printf("%c", a.stringTerm[i]);
    //     }
    //     printf("\n");
    //     printf("A decimal =");
    //     for (int i = 0; i < a.decimalCount; i++) {
    //         printf("%d", a.decimals[i]);
    //     }
    //     printf("\n");
    //     printf("B string =");
    //     for (int i = 0; i < num_vars; i++) {
    //         printf("%c", b.stringTerm[i]);
    //     }
    //     printf("\n");
    //      printf("B decimal =");
    //     for (int i = 0; i < b.decimalCount; i++) {
    //         printf("%d", b.decimals[i]);
    //     }
    //     printf("\n");
    //     printf("DIFFERENCE = %d", *difference);
    //     printf("\n");
    // }


    return res;
}

__global__ void setBinary(int* givens, group* groupings, int num_vars, int* numRemaining, term* remaining, int num_terms) {
    for (int a = 0; a < num_terms; a += BLOCKSIZE) {
        int index = blockIdx.x * blockDim.x + threadIdx.x + a;

        if (index >= num_terms) {
            return;
        }
        // Lock binLock;
        int num_ones = 0;
        //count the ones in the current minterm
        term init_term;

        char strTerm[MAX_VARS + 1];
        int decTerm = givens[index];

        //assign each term to a grouping (depending on how may 1's there are)
        for (int j = num_vars - 1; j >= 0; j--) {
            if (decTerm & 1) {
                strTerm[j] = '1';
                num_ones++;
            } else {
                strTerm[j] = '0';
            }
            decTerm >>= 1;
        }

        //add to the dictionary
        for (int i = 0; i <= num_vars; i++) {
            init_term.stringTerm[i] = strTerm[i]; // Manual copy instead of strcpy
        }
        
        // printf("DECIMAL TO BE ADDED %d", decTerm);
        init_term.decimals[0] = givens[index];
        init_term.decimalCount = 1;
        init_term.tickmark = 0;
        init_term.groupNum = num_ones;
        
        groupings[num_ones].terms[atomicAdd(&groupings[num_ones].size, 1)] = init_term;
        remaining[atomicAdd(numRemaining, 1)] = init_term;
    }
    return;
}

__global__ void findPrimeImplicants(int num_vars, group* groupings, term* matched, term* primeImp,
                               group* nextGroupings, term* unmatched, term* remaining, term* nextRemaining, 
                               term* primeImplicants, char uniqueStringTerms[MAX_TERMS][MAX_VARS + 1],
                               int* numRemaining, int* nextRemainingId, int* numMatched, int* numUnmatched, 
                               int* numPrimeImp, int* notDone, int num_terms){

    for (int a = 0; a < num_terms; a += BLOCKSIZE) {
        int index = blockIdx.x * blockDim.x + threadIdx.x + a;

        if (index >= num_terms) {
            return; 
        }
        //printf("groupings initialized\n");

            while (true) {
    
                __syncthreads();

                if (index % BLOCKSIZE == 0) {
                    *notDone = 0;
                }

                __syncthreads();

                if (index < *numRemaining) {
                    
                    term curr = remaining[index];
                    int numOnes = curr.groupNum;

                    if (numOnes < num_vars) {
                
                        for (int i = 0; i < groupings[numOnes + 1].size; i++) {
                    
                            int differenceValue = 0;       // Create a local integer
                            int* difference = &differenceValue;  // Point to the local integer

                        
                            term match = setMatch(curr, groupings[numOnes + 1].terms[i], num_vars, difference);
                
                            if (*difference == 1){ //matched
                        
                                atomicAdd(notDone, 1);
                                int onesGroup = match.groupNum;

                                nextRemaining[atomicAdd(nextRemainingId, 1)] = match;
                    
                                nextGroupings[onesGroup].terms[atomicAdd(&nextGroupings[onesGroup].size, 1)] = match;

                                matched[atomicAdd(numUnmatched, 1)] = match;

                                //record that terms have been matched for when they are compared again
                                curr.tickmark = 1;
                                groupings[numOnes+1].terms[i].tickmark = 1;
            
                            } 
                            
                            else { //not matched
                                //if we are comparing to the last group and one is left unmatched
                                if (numOnes == num_vars - 1 && groupings[numOnes+1].terms[i].tickmark == 0) {
                                    unmatched[atomicAdd(numUnmatched, 1)] = groupings[numOnes+1].terms[i];
                                } 
                            }
                        }
                    }

                    if (curr.tickmark == 0) {
                        unmatched[atomicAdd(numUnmatched, 1)] = curr;
                    }
                }

                if (!(*notDone)) {
                    break;
                }
                

                if (index < *numMatched) {
                    primeImp[index] = matched[index];
                    *numPrimeImp = *numMatched;
                    *numMatched = 0;
                }

                //remake groupings for the next iteration
                if (index < num_vars + 1) {
                    groupings[index] = nextGroupings[index];
                    nextGroupings[index].size = 0;
                }
                if (index < *nextRemainingId) {
                    remaining[index] = nextRemaining[index];
                    *numRemaining = *nextRemainingId;
                    *nextRemainingId = 0;
                }

            }   
        
        //the groupings should now hold the minimal matchings 
        
        if (index < *numPrimeImp) {
            term match = primeImp[index];
            primeImplicants[index] = match;
        }

        if (index >= *numPrimeImp && index < *numUnmatched + *numPrimeImp) {
            term temp = unmatched[index];
            primeImplicants[index] = temp;
        }

        if (index == 0) {
            *numPrimeImp = *numMatched + *numUnmatched;
        }

    //      if (index == 0) {
    //     for (int j = 0; j < *numUnmatched + *numPrimeImp; j++) {

    //             printf("  Term %d: %s (Decimals: ", 
    //                    j, 
    //                    primeImplicants[j].stringTerm);
                
    //             // Loop through all decimals for the term
    //             for (int d = 0; d < num_terms; d++) {
    //                 printf("%d", primeImplicants[j].decimals[d]);
    //                 if (d < num_terms - 1) {
    //                     printf(", "); // Add a comma between decimals
    //                 }
    //             }
    //             printf(")\n"); // Close the decimal list
    //     }
    // }
        
        //printf("PRIMEIMPLICANT %d STR IS %s", index, primeImplicants[index].stringTerm);
    }
  return;
}

//Function to find essential prime implicants

__global__ void findEssentialImplicants(int num_vars, term* primeImp, int* numPrimeImp, term* essentialImp, int* numEssentials) {

    for (int a = 0; a < *numPrimeImp; a += BLOCKSIZE) {
        int index = blockIdx.x * blockDim.x + threadIdx.x + a;
        if (index < *numPrimeImp) {
            bool isEssential = false;
            term currImp = primeImp[index];
            for (int j = 0; j < currImp.decimalCount; j++) { // Iterate over all minterms covered by primeImp[i]
                int currentDecimal = currImp.decimals[j];
                bool unique = true;

                // Check if this minterm is covered by any other prime implicant
                for (int k = 0; k < *numPrimeImp; k++) {
                    if (k != index) { // Don't compare the implicant with itself
                        for (int l = 0; l < primeImp[k].decimalCount; l++) {
                            if (currentDecimal == primeImp[k].decimals[l]) {
                                unique = false;
                                break;
                            }
                        }
                    }
                    if (!unique) break; // Exit early if duplicate found
                }

                // If a unique minterm is found, mark the implicant as essential
                if (unique) {
                    isEssential = true;
                    break;
                }
            }

            // Add the implicant to the list of essential implicants if marked
            if (isEssential) {
                essentialImp[atomicAdd(numEssentials, 1)] = currImp;
            }
        }
    
    }

    return;
}

term* launchPrimeImplicants(int num_vars, int num_terms, int* givens) {
    // Allocate memory for host variables
    group* groupings = (group*)malloc(sizeof(group) * (num_vars + 1));
    group* nextGroupings = (group*)malloc(sizeof(group) * (num_vars + 1));
    term* primeImplicants = (term*)malloc(sizeof(term) * num_terms);

    if (!groupings || !nextGroupings || !primeImplicants) {
        std::cerr << "Memory allocation failed on host." << std::endl;
        return nullptr;
    }

    // Initialize group sizes
    for (int i = 0; i <= num_vars; i++) {
        groupings[i].size = 0;
        nextGroupings[i].size = 0;
    }

    // Allocate device memory
    int *d_givens;
    group *d_groupings, *d_nextGroupings;
    term *d_matched, *d_primeImp, *d_unmatched, *d_remaining, *d_nextRemaining, *d_primeImplicants, *d_essentialImp;
    char (*d_uniqueStringTerms)[MAX_VARS + 1];
    int *d_numRemaining, *d_nextRemainingId, *d_numMatched, *d_numUnmatched, *d_numPrimeImp, *d_notDone, *d_numEssentials;

    cudaMalloc(&d_givens, sizeof(int) * num_terms);
    cudaMalloc(&d_groupings, sizeof(group) * (num_vars + 1));
    cudaMalloc(&d_nextGroupings, sizeof(group) * (num_vars + 1));
    cudaMalloc(&d_matched, sizeof(term) * num_terms);
    cudaMalloc(&d_primeImp, sizeof(term) * num_terms);
    cudaMalloc(&d_unmatched, sizeof(term) * num_terms);
    cudaMalloc(&d_remaining, sizeof(term) * num_terms);
    cudaMalloc(&d_nextRemaining, sizeof(term) * num_terms);
    cudaMalloc(&d_primeImplicants, sizeof(term) * num_terms);
    cudaMalloc(&d_uniqueStringTerms, sizeof(char) * MAX_TERMS * (MAX_VARS + 1));
    cudaMalloc(&d_numRemaining, sizeof(int));
    cudaMalloc(&d_nextRemainingId, sizeof(int));
    cudaMalloc(&d_numMatched, sizeof(int));
    cudaMalloc(&d_numUnmatched, sizeof(int));
    cudaMalloc(&d_numPrimeImp, sizeof(int));
    cudaMalloc(&d_notDone, sizeof(int));
    cudaMalloc(&d_numEssentials, sizeof(int));
    

    // Check CUDA memory allocation success
    if (cudaPeekAtLastError() != cudaSuccess) {
        std::cerr << "CUDA memory allocation failed." << std::endl;
        return nullptr;
    }

    // Initialize device variables
    int initial_notDone = 0;
    int initial_numRemaining = 0;
    int initial_numMatched = 0;
    int initial_numUnmatched = 0;
    int initial_numPrimeImp = 0;
    int initial_numEssential = 0;

    cudaMemcpy(d_givens, givens, sizeof(int) * num_terms, cudaMemcpyHostToDevice);
    cudaMemcpy(d_groupings, groupings, sizeof(group) * (num_vars + 1), cudaMemcpyHostToDevice);
    cudaMemcpy(d_numRemaining, &initial_numRemaining, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_numMatched, &initial_numMatched, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_numUnmatched, &initial_numUnmatched, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_numPrimeImp, &initial_numPrimeImp, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_notDone, &initial_notDone, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(d_numEssentials, &initial_numEssential, sizeof(int), cudaMemcpyHostToDevice);

    // Configure kernel launch parameters
    int threads_per_block = 256;
    int num_blocks = (num_terms + threads_per_block - 1) / threads_per_block;

    setBinary<<<32,32>>>(d_givens, d_groupings, num_vars, d_numRemaining, d_remaining, num_terms);

    // Launch the CUDA kernel
    findPrimeImplicants<<<32,32>>>(
        num_vars, d_groupings, d_matched, d_primeImp, d_nextGroupings, d_unmatched,
        d_remaining, d_nextRemaining, d_primeImplicants, d_uniqueStringTerms,
        d_numRemaining, d_nextRemainingId, d_numMatched, d_numUnmatched, d_numPrimeImp,
        d_notDone, num_terms);

    printf("Found Prime Implicants\n");


    int numPrimeImp;
    cudaMemcpy(&numPrimeImp, d_numPrimeImp, sizeof(int), cudaMemcpyDeviceToHost);
    term* essentialImp = (term*)malloc(sizeof(term) * numPrimeImp);
    cudaMalloc(&d_essentialImp, sizeof(term) * numPrimeImp);
    cudaMemcpy(d_essentialImp, essentialImp, sizeof(term) * numPrimeImp, cudaMemcpyHostToDevice);

    findEssentialImplicants<<<32,32>>>(
        num_vars, d_primeImplicants, d_numPrimeImp, d_essentialImp, d_numEssentials);

    printf("Found Essential Prime Implicants");

    // Check for kernel launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        std::cerr << "CUDA kernel launch error: " << cudaGetErrorString(err) << std::endl;
        return nullptr;
    }

    // Synchronize device
    cudaDeviceSynchronize();

    // Copy results back to host
    cudaMemcpy(essentialImp, d_essentialImp, sizeof(term) * numPrimeImp, cudaMemcpyDeviceToHost);

    // Free device memory
    cudaFree(d_givens);
    cudaFree(d_groupings);
    cudaFree(d_nextGroupings);
    cudaFree(d_matched);
    cudaFree(d_primeImp);
    cudaFree(d_unmatched);
    cudaFree(d_remaining);
    cudaFree(d_nextRemaining);
    cudaFree(d_primeImplicants);
    cudaFree(d_uniqueStringTerms);
    cudaFree(d_numRemaining);
    cudaFree(d_nextRemainingId);
    cudaFree(d_numMatched);
    cudaFree(d_numUnmatched);
    cudaFree(d_numPrimeImp);
    cudaFree(d_notDone);
    cudaFree(d_essentialImp);
    cudaFree(d_numEssentials);

    // Free host memory
    free(groupings);
    free(nextGroupings);
    free(primeImplicants);

    return essentialImp;
}

term* seqQuineMcCluskey(int argc, char *argv[], int num_vars, int* givens, int num_terms){
    term* primeImps = launchPrimeImplicants(num_vars, num_terms, givens);
    return primeImps;
}


int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  //parse the command line to get the number of variables and all the ones terms
  std::vector<int> given;
  int num_vars = parseInput(argc, argv, given);
  int num_terms = given.size();
  int* givens = (int*)malloc(sizeof(int) * num_terms);

  for (int i = 0; i < num_terms; i++) {
    givens[i] = given[i];
  }
  printf("PARSED\n");

  term* primeImp = seqQuineMcCluskey(argc, argv, num_vars, givens, num_terms);

//   for (int i = 0; i < primeImp.size(); i++) {
//         std::cout << "Implicant " << i + 1 << ": ";
//         for (int dec : primeImp[i].decimals) {
//             std::cout << dec << " ";
//         }
//         std::cout << '\n';
//   }
    
  const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
  std::cout << "Computation time (sec): " << compute_time << '\n';

}
