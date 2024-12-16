/**
 * Sequential K-Map Reduction via QuineMcCluskey algorithm
 * Emily Allendorf(eallendo), Andrew Wang(herongw)
 */

#include <algorithm>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <chrono>
#include <string>
#include <cstring>
#include <vector>

#include <unistd.h>
#include <omp.h>
#include <math.h>

// struct group {
//     std::vector<std::string> terms;
//     int size;
// };

struct term {
    std::string stringTerm;
    std::vector<int> decimals;
    int tickmark;
};

struct group {
    std::vector<term> terms;
    int size; //the number of terms in the group
};

int setMatch(const term& a, const term& b, int num_vars, term& res, int& difference) {
    int ones = 0;
    for (int i = 0; i < num_vars; i++){
        if (a.stringTerm[i] != b.stringTerm[i]) {
            difference++;
            res.stringTerm += "_";
        } else {
            const char bin = a.stringTerm[i];
            // res.stringTerm += a.stringTerm[i];
            // ones += std::stoi(a.stringTerm[i]);
            res.stringTerm += bin;
            if (bin == '1') ones += 1;
        }
    }
    //populate the result decimals with all of the minterms involved
    res.decimals.insert(res.decimals.end(), a.decimals.begin(), a.decimals.end());
    res.decimals.insert(res.decimals.end(), b.decimals.begin(), b.decimals.end());
    res.tickmark = 0;

    // std::cout << "A: " << a.stringTerm << " (Decimals: ";
    // for (int dec : a.decimals) std::cout << dec << " ";
    // std::cout << ")\n";
    // std::cout << "B: " << b.stringTerm << " (Decimals: ";
    // for (int dec : b.decimals) std::cout << dec << " ";
    // std::cout << ")\n";
    //   std::cout << "Resulting Term: " << res.stringTerm << "\n";
    // std::cout << "Difference: " << difference << "\n";

    return ones;
}


std::vector<term> findPrimeImplicants(int num_vars, std::vector<int>givens){
    //initialize the groupings vector
    std::vector<group> groupings(num_vars+1);
    //set the initial groups to size 0
    for (int i = 0; i < num_vars+1; i++){
        groupings[i].size = 0;
    }
    int num_terms = givens.size();
    //assign each term to a grouping (depending on how may 1's there are)
    for (int i = 0; i < num_terms; i++) {
        term init_term;
        //count the ones in the current minterm
        int num_ones = 0;
        std::string strTerm = "";
        int decTerm = givens[i];
        for (int j = 0; j < num_vars; j++) {
            if (givens[i] & 1){
                // strTerm += "1";
                strTerm.insert(0, "1");
                num_ones = num_ones + (givens[i] & 1);
            } else {
                // strTerm += "0";
                strTerm.insert(0, "0");
            }
            givens[i] >>= 1;
        }
    
        //add to the dictionary
        init_term.stringTerm = strTerm;
        // printf("DECIMAL TO BE ADDED %d", decTerm);
        init_term.tickmark = 0;
        init_term.decimals.push_back(decTerm);
        groupings[num_ones].terms.push_back(init_term);
        groupings[num_ones].size++;
    }
    
    // printf("groupings initialized\n");

    //vector to keep track of the unmatched terms
    std::vector<term> unmatched;

    bool doneMatching;
    //TODO: check the number of iterations there should be
    for (int iters = 0; iters < num_vars-1; iters++){
        
        //build a vector for the next groupings
        std::vector<group> nextGroupings(num_vars+1);
        for (int i = 0; i < num_vars+1; i++){
            nextGroupings[i].size = 0;
        }
        std::vector<term> nextUnmatched;
        
        doneMatching = true;
        //compare groupings num_vars times (there are num_vars+1 groups)
        for (int i = 0; i < num_vars; i++) {
            //compare adjacent groups
            for (int j = 0; j < groupings[i].size; j++) {
                for (int k = 0; k < groupings[i+1].size; k++) {
                    int difference = 0;
                    term match; 
                    int onesGroup = setMatch(groupings[i].terms[j], groupings[i+1].terms[k], num_vars, match, difference);
                    if (difference == 1){ //matched
                        doneMatching = false;
                        nextGroupings[onesGroup].size++;
                        nextGroupings[onesGroup].terms.push_back(match);
                        //record that terms have been matched for when they are compared again
                        groupings[i].terms[j].tickmark = 1;
                        groupings[i+1].terms[k].tickmark = 1;
                    } else { //not matched
                        //if we are comparing to the last group and one is left unmatched
                        if (i == num_vars - 1 && groupings[i+1].terms[k].tickmark == 0) {
                            nextUnmatched.push_back(groupings[i+1].terms[k]);
                        } 
                    }
                }
                if (groupings[i].terms[j].tickmark == 0) {
                    nextUnmatched.push_back(groupings[i].terms[j]);
                }
            }
        }   
        //we are done if there are no new match
        if (doneMatching) break;
        //remake groupings for the next iteration
        groupings = nextGroupings;
        unmatched = nextUnmatched;
    
    }
    
    //the groupings should now hold the minimal matchings 
    std::vector<term> primeImplicants;
    std::vector<std::string> uniqueStringTerms;
    for (int i = 0; i < num_vars+1; i++) {
        for (int j = 0; j < groupings[i].size; j++) {
            //check if the string term is unique to prevent duplicates
            if (!(std::find(uniqueStringTerms.begin(), uniqueStringTerms.end(), 
                groupings[i].terms[j].stringTerm) != uniqueStringTerms.end())) {
                    primeImplicants.push_back(groupings[i].terms[j]);
                    uniqueStringTerms.push_back(groupings[i].terms[j].stringTerm);
            }
        }
    }
    
    std::cout << "Prime Implicants matched:" << std::endl;
    for (const term& implicant : primeImplicants) {
        std::cout << "String: " << implicant.stringTerm << ", Decimals: ";
        for (int dec : implicant.decimals) {
            std::cout << dec << " ";
        }
        std::cout << std::endl;
    }

    //add in the unmatched ones from the previous steps
    int numUnmatched = unmatched.size();
    for (int i = 0; i < numUnmatched; i++){
        if (!(std::find(uniqueStringTerms.begin(), uniqueStringTerms.end(), 
            unmatched[i].stringTerm) != uniqueStringTerms.end())) {
                    primeImplicants.push_back(unmatched[i]);
                    uniqueStringTerms.push_back(unmatched[i].stringTerm);
                } 
    }

//     std::cout << "Prime Implicants:" << std::endl;
// for (const term& implicant : primeImplicants) {
//     std::cout << "String: " << implicant.stringTerm << ", Decimals: ";
//     for (int dec : implicant.decimals) {
//         std::cout << dec << " ";
//     }
//     std::cout << std::endl;
// }
  
  return primeImplicants;
}


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
    // printf("GOTHERE1111");
//     std::cout << "Number of variables: " << num_vars << '\n';
// std::cout << "Minterms: ";
// for (int term : givens) {
//     std::cout << term << " ";
// }
// std::cout << '\n';

    return num_vars;
}

// void findEssenImplicants(int num_vars, int num_terms, std::vector<term>primeImp){
//     std::vector<term> essentialImp; 
//     int numPrime = primeImp.size();
//     for (int i = 0; i < numPrime; i++) { 
//         int numDec = primeImp[i].decimals.size();
//         for (int j = 0; j < numDec, j++) {//ASSUMING THERE IS ALWAYS AT LEAST ONE PRIME IMPLICANT
//             bool unique = true;
//             for (int k = i+1; k < numPrime; k++) {
//                 if (primeImp[i].decimal[j] == primeImp[k].decimal[j]) {
//                     unique = false;
//                 }
//             }
//             if (unique) {
//                 essentialImp.terms.push_back(primeImp[i]);
//             }
//         }
//     }
    
//     return essentialImp;
// }

// Function to find essential prime implicants
std::vector<term> findEssentialImplicants(int num_vars, const std::vector<term>& primeImp) {
    std::vector<term> essentialImp; 
    int numPrime = primeImp.size();

    for (int i = 0; i < numPrime; i++) { 
        bool isEssential = false;
        for (int j = 0; j < primeImp[i].decimals.size(); j++) { // Iterate over all minterms covered by primeImp[i]
            int currentDecimal = primeImp[i].decimals[j];
            bool unique = true;

            // Check if this minterm is covered by any other prime implicant
            for (int k = 0; k < numPrime; k++) {
                if (k != i) { // Don't compare the implicant with itself
                    for (int l = 0; l < primeImp[k].decimals.size(); l++) {
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
            essentialImp.push_back(primeImp[i]);
        }
    }
       std::cout << "ESSENTIAL Implicants:" << std::endl;
for (const term& implicant : essentialImp) {
    std::cout << "String: " << implicant.stringTerm << ", Decimals: ";
    for (int dec : implicant.decimals) {
        std::cout << dec << " ";
    }
    std::cout << std::endl;
}

    return essentialImp; // Return the vector of essential implicants
}

std::vector<term> seqQuineMcCluskey(int argc, char *argv[], int num_vars, std::vector<int>& givens){
    std::vector<term> primeImps = findPrimeImplicants(num_vars, givens);
    printf("found prime implicants\n");
    findEssentialImplicants(num_vars, primeImps);
    return primeImps;

}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();

  const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
  std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';

  const auto compute_start = std::chrono::steady_clock::now();

  //parse the command line to get the number of variables and all the ones terms
  std::vector<int> givens;
  int num_vars = parseInput(argc, argv, givens);
  printf("PARSED\n");
  
  std::vector<term> primeImp = seqQuineMcCluskey(argc, argv, num_vars, givens);

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