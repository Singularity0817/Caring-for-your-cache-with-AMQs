#include <iostream>
#include <fstream>
#include <cstring>
#include <string>
#include <cstdlib>

#include "blocked_bloom_filter.h"



int main(int argc, char **argv)
{
    if(argc < 2)
    {
        std::cerr << "Invalid command.\n";
        exit(1);
    }

    if(!strcmp(argv[1], "build"))
    {
        std::string keyFile;
        double fpr;
        uint64_t n;
        std::string outputFile;

        
        for(int i = 2; i < 10; i += 2)
            if(!strcmp(argv[i], "-k"))
                keyFile = std::string(argv[i + 1]);
            else if(!strcmp(argv[i], "-f"))
                fpr = atof(argv[i + 1]);
            else if(!strcmp(argv[i], "-n"))
                n = atoll(argv[i + 1]);
            else if(!strcmp(argv[i], "-o"))
                outputFile = std::string(argv[i + 1]);
            else
            {
                std::cerr << "Invalid parameter.\n";
                exit(1);
            }


        blocked_bloom_filter bf(keyFile, n, fpr, outputFile);
    }
    else if(!strcmp(argv[1], "query"))
    {
        std::string bbfFile;
        std::string queryFile;

        for(int i = 2; i < 6; i += 2)
            if(!strcmp(argv[i], "-i"))
                bbfFile = std::string(argv[i + 1]);
            else if(!strcmp(argv[i], "-q"))
                queryFile = std::string(argv[i + 1]);
            else
            {
                std::cerr << "Invalid parameter.\n";
                exit(1);
            }        

        
        blocked_bloom_filter bbf(bbfFile);

        std::vector<bool> result;
        bbf.query(queryFile, result);

        for(uint64_t i = 0; i < result.size(); ++i)
            std::cout << "query" << (i + 1) << ":" << (result[i] ? "Y" : "N") << "\n";
    }
    else
    {
        std::cerr << "Invalid command.\n";
    }
    

    return 0;
}