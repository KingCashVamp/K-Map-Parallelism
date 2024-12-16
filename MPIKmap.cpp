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

#include <mpi.h>
#include <Kmap.h>

// #define ROOT 0

//a ring message passing helper to send one integer (aka a size)
// std::vector<std::vector<term>> data_ring(std::vector<int> headerCounts, std::vector<term> data, int size_mult, std::vector<group> groupings, const int nproc, int pid, int nextid, int previd){
//     // std::vector<term> allTerms(total_minterms);
//     std::vector<std::vector<term>> termBuffers(nproc);
//     printf("counts[%d] = %d, div by mult = %d\n", pid, headerCounts[pid], headerCounts[pid]/size_mult);
//     for (int i = 0; i < nproc; i++) {
//         termBuffers[i].resize(headerCounts[i]/size_mult);
//     }
    
//     termBuffers[pid] = data;

//     printf("size of data: %ld, length: %ld, size of data[0]: %ld\n", sizeof(data), data.size(), sizeof(data[0]));

//     MPI_Request sendData;
//         //      24       256                1         
//     MPI_Isend(data.data(), headerCounts[pid], MPI_CHAR, nextid, 0, MPI_COMM_WORLD, &sendData);

//     for (int i = 0; i < nproc - 1; i++) {
//       int recvid = previd - i;
//       if (recvid < 0) {
//         recvid += nproc;
//       }
    
//     //   if (pid == 1) {
//         MPI_Recv(termBuffers[recvid].data(), headerCounts[recvid], MPI_CHAR, previd, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//       printf("%d:  p %d recived from %d into slot %d\n", headerCounts[recvid], pid, previd, recvid);
//     //   }

//       MPI_Request propData;
//       if (nextid != pid && nproc > 2 && (i != nproc - 1)) {
//         MPI_Isend(termBuffers[recvid].data(), headerCounts[recvid], MPI_CHAR, nextid, 0, MPI_COMM_WORLD, &propData);
//       }
//     }

//     //should we update the groupings in here?
    
//     return termBuffers;
// }

std::vector<std::vector<partialTerm>> data_ring(std::vector<long> headerCounts, std::vector<partialTerm> data, int size_mult, std::vector<group> groupings, const int nproc, int pid, int nextid, int previd, MPI_Datatype MPI_PTERM){
    // std::vector<term> allTerms(total_minterms);
    std::vector<std::vector<partialTerm>> termBuffers(nproc);
    // termBuffers[pid].resize(headerCounts[pid]/size_mult);
    // termBuffers[pid] = data;

    
    printf("counts[%d] = %ld, div by mult = %ld\n", pid, headerCounts[pid], headerCounts[pid]/size_mult);
    for (int i = 0; i < nproc; i++) {
        termBuffers[i].resize(headerCounts[i]/size_mult);
        if (i == pid) termBuffers[i] = data;
        long num_terms = headerCounts[i]/size_mult;
        for (int j = 0; j < num_terms; j++){
            printf("p%d[%d]: %s\n", pid, i, termBuffers[i][j].stringTerm);
        }
        // if (headerCounts[i]%size_mult != 0) printf("asldJKSHGULIUR\n");
    }

    std::vector<partialTerm> sendData = termBuffers[pid];
    std::vector<partialTerm>  recvData;

    for (int i = 0; i < nproc-1; i++) {
      int recvid = (previd - i + nproc) % nproc;
    //   if (recvid < 0) {
    //     recvid += nproc;
    //   }
      int sendid = (pid + i) % nproc;
    //   if (sendid > nproc) {
    //     sendid -= nproc;
    //   }

      MPI_Request propData;
    //   if (nproc > 2 && (i != nproc - 1)) {
        printf("i = %d; p%d sending out slot %d's data to p%d\n",i,  pid, sendid,nextid);

        // MPI_Isend(termBuffers[sendid].data(), headerCounts[sendid], MPI_CHAR, nextid, 0, MPI_COMM_WORLD, &propData);
        MPI_Isend(sendData.data(), sendData.size(), MPI_PTERM, nextid, sendid, MPI_COMM_WORLD, &propData);
        printf("p%d sent out slot %d's data to p%d\n", pid, sendid,nextid);
    //   }
    
      printf("i = %d; p %d reciving from %d into slot %d\n  recvid size: %ld, message size: %ld\n", i, pid, previd, recvid, headerCounts[recvid], headerCounts[previd]);
            
      recvData.resize(headerCounts[recvid]/size_mult);
      printf("***: %ld\n", headerCounts[recvid]/size_mult);

      MPI_Request recvReq;
      MPI_Irecv(recvData.data(), recvData.size(), MPI_PTERM, previd, recvid, MPI_COMM_WORLD, &recvReq);
    //   MPI_Wait(&recvReq, MPI_STATUS_IGNORE);
    //   MPI_Recv(termBuffers[recvid].data(), headerCounts[recvid], MPI_CHAR, previd, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      printf("p %d recived from %d into slot %d\n  recvid size: %ld, message size: %ld\n", pid, previd, recvid, headerCounts[recvid], headerCounts[previd]);
    
    sendData = recvData;
    for (int k = 0; k < recvData.size(); k++){
        printf("recvid: %d, [%d] = %s\n", recvid, k, recvData[k].stringTerm);
    }
    termBuffers[recvid] = recvData;

    }
    
    // termBuffers[pid] = data;
    // printf("init p%d: ")

    // printf("size of data: %ld, length: %ld, size of data[0]: %ld\n", sizeof(data), data.size(), sizeof(data[0]));

    // printf("size: %ld\n", headerCounts[pid]);
    // MPI_Request sendData;
    //     //      24       256                1         
    // if (nextid != pid) {
    //     printf("p%d sent out slot %d's data to p%d\n", pid, pid,nextid);
    //     MPI_Isend(data.data(), headerCounts[pid], MPI_CHAR, nextid, pid, MPI_COMM_WORLD, &sendData);
    // }

    // for (int i = 0; i < nproc - 1; i++) {
    //   int recvid = previd - i;
    //   if (recvid < 0) {
    //     recvid += nproc;
    //   }
    
    // //   if (pid == 1) {

    //     MPI_Recv(termBuffers[recvid].data(), headerCounts[recvid], MPI_CHAR, previd, recvid, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //   printf("p %d recived from %d into slot %d\n  recvid size: %ld, message size: %ld\n", pid, previd, recvid, headerCounts[recvid], headerCounts[previd]);
    // //   }

    //   MPI_Request propData;
    //   if (nextid != pid && nproc > 2 && (i != nproc - 1)) {
    //     MPI_Isend(termBuffers[recvid].data(), headerCounts[recvid], MPI_CHAR, nextid, recvid, MPI_COMM_WORLD, &propData);
    //     printf("p%d sent out slot %d's data to p%d\n", pid, recvid,nextid);
    //   }
    // }

    //should we update the groupings in here?
    
    return termBuffers;
}

//a ring message passing helper to send one integer (aka a size)
std::vector<long> header_ring(long header, const int nproc, int pid, int nextid, int previd){
    //send a header that has the assignedTerms size to all the other processors
    MPI_Request sendHeader;
    MPI_Isend(&header, 1, MPI_LONG, nextid, 0, MPI_COMM_WORLD, &sendHeader);
    //recieve the header information in a ring
    std::vector<long> allHeaders(nproc);
    allHeaders[pid] = header;
    for (int i = 0; i < nproc - 1; i++) {
      int recvid = previd - i;
      if (recvid < 0) {
        recvid += nproc;
      }
      //recieve the previous header into the correct place in the headers vector
      MPI_Recv(&allHeaders[recvid], 1, MPI_LONG, previd, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      //propagate the header information
      MPI_Request propHeader;
      if (nextid != pid && nproc > 2 && (i != nproc - 1)) {
        MPI_Isend(&allHeaders[recvid], 1, MPI_LONG, nextid, 0, MPI_COMM_WORLD, &propHeader);
      }
    }

    return allHeaders;
}


//may need to refactor this later
int setMatch(const term& a, const term& b, int num_vars, term& res, int& difference){//, std::vector<int> decimals) {
    int ones = 0;
    for (int i = 0; i < num_vars; i++){
        if (a.stringTerm[i] != b.stringTerm[i]) {
            difference++;
            // res.stringTerm += "_";
            res.stringTerm[i] = '_';
        } else {
            const char bin = a.stringTerm[i];
            // res.stringTerm += a.stringTerm[i];
            // ones += std::stoi(a.stringTerm[i]);
            // res.stringTerm += bin;
            res.stringTerm[i] = bin;
            if (bin == '1') ones += 1;
        }
    }
    //populate the result decimals with all of the minterms involved
    res.stringTerm[num_vars] = '\0';

    res.decimals.insert(res.decimals.end(), a.decimals.begin(), a.decimals.end());
    res.decimals.insert(res.decimals.end(), b.decimals.begin(), b.decimals.end());
    res.tickmark = 0;
    res.groupid = ones;

    // decimals.insert(decimals.end(), a.decimals.begin(), a.decimals.end());
    // decimals.insert(decimals.end(), b.decimals.begin(), b.decimals.end());


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

//let the givens passed in be the assigned minterms for the processor
std::vector<term> findPrimeImplicants(int num_vars, std::vector<int>givens, int total_minterms, int local_minterms, const int nproc, int pid, std::vector<int> counts, MPI_Datatype MPI_PTERM){//, MPI_Datatype MPI_term){
    
    // struct term {
    //     int tickmark;
    //     int groupid;
    //     char[num_vars] stringTerm; //max length of num_vars easy fix the length should be the same size always
    //     // char *stringTerm;
    //     int[1] decimals; //max length would be the number of minterms which is like to long
    
    // };
    
    //logic necessary for message passing
    // int nextid, previd;
    // if (pid == nproc - 1) {
    //   nextid = 0;
    //   previd = pid - 1;
    // } 
    // else if (pid == 0) {
    //   nextid = pid + 1;
    //   previd = nproc - 1;
    // } 
    // else {
    //   nextid = pid + 1;
    //   previd = pid - 1;
    // }  

    int nextid, previd;

// For the first process (pid == 0)
if (pid == 0) {
    nextid = 1;                // Process 0 sends to process 1
    previd = nproc - 1;        // Process 0 receives from the last process (nproc - 1)
}
// For the last process (pid == nproc - 1)
else if (pid == nproc - 1) {
    nextid = 0;                // The last process sends to process 0 (closing the ring)
    previd = pid - 1;          // The last process receives from the second-to-last process
}
// For all other intermediate processes
else {
    nextid = pid + 1;          // Send to the next process
    previd = pid - 1;          // Receive from the previous process
}
    
    
    //initialize the groupings vector
    std::vector<group> groupings(num_vars+1);
    //initialize a local vector of terms the processor will be responsible for
    std::vector<partialTerm> assignedTerms;
    
    int terms_count = 0;
    // printf("successful, pid: %d, header: %d\n", pid, header);
    //set the initial groups to size 0
    for (int i = 0; i < num_vars+1; i++){
        groupings[i].size = 0;
    }
    
    //make a vector of terms for the processor
    int num_terms = local_minterms;
    //assign each term to a grouping (depending on how may 1's there are)
    std::vector<int> assignedDecimals;
    for (int i = 0; i < num_terms; i++) {
        partialTerm init_term;
        //count the ones in the current minterm
        int num_ones = 0;
        // std::string strTerm = "";
        int dec = givens[i];
        for (int j = 0; j < num_vars; j++) {
            if (givens[i] & 1){
                // strTerm.insert(0, "1");
                init_term.stringTerm[num_vars-1-j] = '1';
                num_ones = num_ones + (givens[i] & 1);
            } else {
                // strTerm.insert(0, "0");
                init_term.stringTerm[num_vars-1-j] = '0';
            }
            givens[i] >>= 1;
        }
        init_term.stringTerm[num_vars] = '\0'; //is this not implicit?
    
        //add to the local vector of terms
        // init_term.stringTerm = strTerm;
        assignedDecimals.push_back(dec);
        // decTerms.push_back(dec);
        // decimals.push_back(decTerms);
        init_term.tickmark = 0;
        init_term.groupid = num_ones;
        assignedTerms.push_back(init_term);
        terms_count++;
        // groupings[num_ones].terms.push_back(init_term);
        // groupings[num_ones].size++;
    }
//delimiter?
    // int size_mult = ((sizeof(char)*(num_vars+1)) +
    //                  (sizeof(int)*2        ) +
    //                  (sizeof(int)*(1)      ));
    // printf("0: %ld\n", sizeof(assignedTerms[0]));
    // printf("1: %ld\n", sizeof(assignedTerms[1]));
    // printf("2: %ld\n", sizeof(assignedTerms[2]));
    // printf("3: %ld\n", sizeof(assignedTerms[3]));
    long size_mult = sizeof(assignedTerms[0]); //should all be contiguous memory
    
    // int getBackSize = allSizes[pid]/size_mult;
    // printf("getBackSize: %d\n", getBackSize);
    // for (int i = 0; i < nproc; i++){
        // for (int j = 0; j < getBackSize; j++) {
    //     for (int j = 0; j < local_minterms; j++) {
    //         // if (pid == 0) printf("pid %d: [%d][%d], %s\n", pid, i ,j, allData[i][j].stringTerm);
    //         printf("pid %d: [%d][%d], %s\n", pid, pid ,j, assignedTerms[j].stringTerm);
    //         // printf("%d\n", allData[i][j].decimals[0]);
    //         // groupings[allData[i][j].groupid].terms.push_back(allData[i][j]);
    //     }
    // // }

    long local_size = (long)terms_count*size_mult;

    //this call is lowkey unnescessary comm but good for testing
    std::vector<long> allSizes =  header_ring(local_size, nproc, pid, nextid, previd);
    printf("pid %d: header: %ld \n", pid, allSizes[pid]);

    

    std::vector<std::vector<partialTerm>> allData = data_ring(allSizes, assignedTerms, size_mult, groupings, nproc, pid, nextid, previd, MPI_PTERM);

    //decimals_ring
    MPI_Barrier(MPI_COMM_WORLD);
    // printf("---> %d\n", allData[0][0].decimals[0]);
    
    int getBackSize = allSizes[pid]/size_mult;
    // printf("getBackSize: %d\n", getBackSize);
    for (int i = 0; i < nproc; i++){
        for (int j = 0; j < getBackSize; j++) {
    //     for (int j = 0; j < local_minterms; j++) {
            if (pid == 0) printf("pid %d: [%d][%d], %s\n", pid, i ,j, allData[i][j].stringTerm);
    //         printf("pid %d: [%d][%d], %s\n", pid, pid ,j, assignedTerms[j].stringTerm);
    //         // printf("%d\n", allData[i][j].decimals[0]);
    //         // groupings[allData[i][j].groupid].terms.push_back(allData[i][j]);
        }
    }


    // for (int gr = 0; gr < groupings.size(); gr++){
    //     for (int t = 0; t < groupings[gr].terms.size(); t++){
    //         printf("pid: %d, %s %d\n", pid, groupings[gr].terms[t].stringTerm.c_str(), groupings[gr].terms[t].decimals[0]);
    //     }
    // }

    // for (int i = 0; i < nproc; i++){
    //     printf("p %d: headers[%d] = %d\n", pid, i, allSizes[i]);
    // }
    // printf("headers[%d] = %d, counts[%d] = %d\n", pid, allHeaders[pid], pid, counts[pid]);
    //end useless ring section

    
    
    
    
    
    
    
    
    
    
    std::vector<term> test;
    return test;

    // printf("groupings initialized\n");

    
    
    
    
    
    
    
    
    
    
    //Sequential code from here on down should work if the above code updates the groupings correctly
    //vector to keep track of the unmatched terms
    
    std::vector<term> unmatched;

    bool doneMatching = false;
    //TODO: check the number of iterations there should be
    // for (int iters = 0; iters < num_vars-1; iters++){
    while (!doneMatching) {
        
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
                        // printf("%s\n", groupings[i].terms[j].stringTerm.c_str());
                        doneMatching = false;
                        nextGroupings[onesGroup].size++;
                        nextGroupings[onesGroup].terms.push_back(match);
                        //record that terms have been matched for when they are compared again
                        groupings[i].terms[j].tickmark = 1;
                        groupings[i+1].terms[k].tickmark = 1;
                    } else { //not matched
                        //if we are comparing to the last group and one is left unmatched
                        if (i == num_vars - 1 && (j == groupings[i].size - 1) && groupings[i+1].terms[k].tickmark == 0) {
                            nextUnmatched.push_back(groupings[i+1].terms[k]);
                        } 
                    }
                }
                // printf("%s\n", groupings[i].terms[j].stringTerm.c_str());
                // printf("%d\n", groupings[i].terms[j].tickmark);
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
    
    // std::cout << "Prime Implicants matched:" << std::endl;
    // for (const term& implicant : primeImplicants) {
    //     std::cout << "String: " << implicant.stringTerm << ", Decimals: ";
    //     for (int dec : implicant.decimals) {
    //         std::cout << dec << " ";
    //     }
    //     std::cout << std::endl;
    // }

    //add in the unmatched ones from the previous steps
    int numUnmatched = unmatched.size();
    // printf("numUnmatched: %d\n", numUnmatched);
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

std::vector<term> seqQuineMcCluskey(int argc, char *argv[], int num_vars, std::vector<int>& givens, int total_minterms, int local_minterms, const int nproc, int pid, std::vector<int> counts, MPI_Datatype MPI_PTERM){
    if (pid == 1) printf("successful, minterms: %d,\n", local_minterms);
    std::vector<term> primeImps = findPrimeImplicants(num_vars, givens, total_minterms, local_minterms, nproc, pid, counts, MPI_PTERM);
    printf("found prime implicants\n");
    findEssentialImplicants(num_vars, primeImps);
    return primeImps;

}

int main(int argc, char *argv[]) {
  const auto init_start = std::chrono::steady_clock::now();
  int pid;
  int nproc;

  // Initialize MPI
  MPI_Init(&argc, &argv);
  // Get process rank
  MPI_Comm_rank(MPI_COMM_WORLD, &pid);
  // Get total number of processes
  MPI_Comm_size(MPI_COMM_WORLD, &nproc);

  std::string input_filename;
  int num_vars;
  int opt;
  while ((opt = getopt(argc, argv, "f:n:")) != -1) {
    switch (opt) {
      case 'f':
        input_filename = optarg;
        break;
      case 'n':
        num_vars = atoi(optarg);
        break;
      default:
        if (pid == ROOT) {
          std::cerr << "Reading opts. Usage: " << argv[0] << " -f input_filename -n num_vars\n";
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
    }
  }
  if (empty(input_filename)) {
    if (pid == ROOT) {
          std::cerr << "Usage: " << argv[0] << " -f input_filename -n num_vars\n";
        }
        MPI_Finalize();
        exit(EXIT_FAILURE);
  }
  
  //make an MPI datatype to match the term struct
  int block_lengths[3] = {1, 1, 5};
  MPI_Datatype block_types[3] = {MPI_INT, MPI_INT, MPI_CHAR};
  MPI_Aint     displacements[3];
  displacements[0] = offsetof(partialTerm, tickmark);
  displacements[1] = offsetof(partialTerm, groupid);
  displacements[2] = offsetof(partialTerm, stringTerm);
  MPI_Datatype MPI_PTERM;

  MPI_Type_create_struct(3, block_lengths, displacements, block_types, &MPI_PTERM);
  MPI_Type_commit(&MPI_PTERM);


  
  //parse the command line to get the number of variables and all the ones terms
  std::vector<int> givens;
  int num_minterms = 0;
  int last = nproc-1;
  if (pid == ROOT) printf("last: %d\n", last);
//   int num_vars = parseInput(argc, argv, givens);
//   printf("PARSED\n");
  //possible optimizatio: parallelize the file reading
  if (pid == ROOT) {
    std::ifstream fin(input_filename);

    if (!fin) {
        std::cerr << "Unable to open file: " << input_filename << ".\n";
        exit(EXIT_FAILURE);
    }
    int dec;
    while (fin >> dec) {
        givens.push_back(dec);
        num_minterms++;
    }
    printf("number of minterms: %d\n", num_minterms);
    // for (int i = 0; i < num_minterms; i++) {
    //     printf("%d\n", givens[i]);
    // }
  }

  MPI_Bcast(&num_minterms, 1, MPI_INT, 0, MPI_COMM_WORLD);

   int minterms_per_proc = num_minterms / nproc;
//    int minterms_last_group = minterms_per_proc + (num_minterms % nproc);
   int leftover =  (num_minterms % nproc);

    std::vector<int> counts(nproc);
    std::vector<int> displs(nproc);
    for (int i = 0; i < nproc; i++) {
        counts[i] = minterms_per_proc + (i < leftover ? 1 : 0);
        displs[i] = (i == 0) ? 0 : displs[i - 1] + counts[i - 1];
    }

    // Allocate local array for each process
    std::vector<int> local_givens(counts[pid]);

    // Scatter the data
    MPI_Scatterv(givens.data(), counts.data(), displs.data(), MPI_INT,
                 local_givens.data(), counts[pid], MPI_INT, 0, MPI_COMM_WORLD);

    // for (int i = 0; i < counts[pid]; i++){
    //     printf("pid: %d, givens[%d]=%d\n", pid, i, local_givens[i]);
    // }
/*  this works at least
//   MPI_Bcast(&num_minterms, 1, MPI_INT, 0, MPI_COMM_WORLD);
//   if (pid != ROOT) {
//     givens.resize(num_minterms);
//   }
// //   printf("first succesful\n");

//   MPI_Bcast(givens.data(), num_minterms, MPI_INT, 0, MPI_COMM_WORLD);
// //   printf("second successful");
*/

  if (pid == ROOT) {
    // for (int i = 0; i < num_minterms; i++) {
    //     printf("%d\n", givens[i]);
    // }
    const double init_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - init_start).count();
    std::cout << "Initialization time (sec): " << std::fixed << std::setprecision(10) << init_time << '\n';
  }

  const auto compute_start = std::chrono::steady_clock::now();
  
//   //for now just have the root processor do everything
//   if (pid == ROOT) {
//     std::vector<term> primeImp = seqQuineMcCluskey(argc, argv, num_vars, givens, num_minterms);
//     // printf("seq successful\n");
//   }


  std::vector<term> primeImp = seqQuineMcCluskey(argc, argv, num_vars, local_givens, num_minterms, counts[pid], nproc, pid, counts, MPI_PTERM);


  if (pid == ROOT) { 
    const double compute_time = std::chrono::duration_cast<std::chrono::duration<double>>(std::chrono::steady_clock::now() - compute_start).count();
    std::cout << "Computation time (sec): " << compute_time << '\n';
  }

  //cleanup
  MPI_Type_free(&MPI_PTERM);
  MPI_Finalize();
}