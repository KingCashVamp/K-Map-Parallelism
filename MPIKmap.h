/**
 * Parallel VLSI Wire Routing via OpenMP
 * Name 1(andrew_id 1), Name 2(andrew_id 2)
 */

// #include <omp.h>
#include <cstdint>
#include <vector>
#include <string>

#define ROOT 0
// #define MAX_DIM;

struct partialTerm {
    int tickmark;
    int groupid;
    char stringTerm[5]; //hardcoded num_vars
    // std::string stringTerm; //max length of num_vars easy fix the length should be the same size always
    // char *stringTerm;
    // std::vector<int> decimals; //max length would be the number of minterms which is like to long
};

// struct fullTerm {
//     term partialTerm;
//     std::vector<int> decimals;
// };

struct term {
    int tickmark;
    int groupid;
    std::string stringTerm; //max length of num_vars easy fix the length should be the same size always
    std::vector<int> decimals; //max length would be the number of minterms which is like to long
};
struct group {
    std::vector<term> terms;
    int size; //the number of terms in the group
};