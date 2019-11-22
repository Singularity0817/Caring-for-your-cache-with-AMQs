#include <chrono>
#include <random>
#include <vector>
#include <unordered_set>

#include "blocked_bloom_filter.h"


void generate_random_list(uint64_t sz, std::vector<uint64_t> &l)
{
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    l.reserve(sz);

    for(uint64_t i = 0; i < sz; ++i)
        l.push_back(std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(rng));
}


void generate_random_list(uint64_t sz, std::vector<uint64_t> &l,
                            std::unordered_multiset<uint64_t> &forbidden)
{
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    for(uint64_t i = 0; i < sz; ++i)
        while(true)
        {
            uint64_t val = std::uniform_int_distribution<uint64_t>(0, std::numeric_limits<uint64_t>::max())(rng);
            if(forbidden.find(val) == forbidden.end())
            {
                l.push_back(val);
                break;
            }
        }
}


void benchmark_bf_all_negative(const uint64_t elemsStart, const uint64_t elemsEnd, const uint64_t elemsStep,
                                const double fprStart, const double fprEnd, const double fprStep,
                                const uint64_t queryCount)
{
    std::vector<uint64_t> keyList;
    std::vector<uint64_t> queryList;

    generate_random_list(queryCount, queryList);
    std::unordered_multiset<uint64_t> querySet(queryList.begin(), queryList.end());

    for(double fpr = fprStart; fpr <= fprEnd; fpr += fprStep)
        for(uint64_t n = elemsStart; n <= elemsEnd; n += elemsStep)
        {
            keyList.clear();
            
            generate_random_list(n, keyList, querySet);
            bloom_filter bf(n, fpr);

            for(auto p = keyList.begin(); p != keyList.end(); ++p)
                bf.insert(*p);

            // bf.dump_metadata();


            double elapsedSecs = 0;
            uint64_t positive = 0;
            
            for(auto p = queryList.begin(); p != queryList.end(); ++p)
            {
                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                positive += bf.query(*p);
                std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
                elapsedSecs += time_span.count();
            }


            std::clog << fpr << "\t" << n << "\t" << elapsedSecs << "\t" << double(positive) / queryCount << "\n";
        }
}



void benchmark_bf_half_negative(const uint64_t elemsStart, const uint64_t elemsEnd, const uint64_t elemsStep,
                                const double fprStart, const double fprEnd, const double fprStep,
                                const uint64_t queryCount)
{
    std::vector<uint64_t> keyList;
    std::vector<uint64_t> queryList;


    for(double fpr = fprStart; fpr <= fprEnd; fpr += fprStep)
        for(uint64_t n = elemsStart; n <= elemsEnd; n += elemsStep)
        {
            keyList.clear();
            queryList.clear();

            generate_random_list(n, keyList);

            std::unordered_multiset<uint64_t> keySet(keyList.begin(), keyList.end());

            for(uint64_t i = 0; i < queryCount / 2; ++i)
                queryList.push_back(keyList[i]);
            
            generate_random_list(queryCount / 2, queryList, keySet);


            bloom_filter bf(n, fpr);

            for(auto p = keyList.begin(); p != keyList.end(); ++p)
                bf.insert(*p);

            // bf.dump_metadata();


            double elapsedSecs = 0;
            uint64_t positive = 0;
            
            for(auto p = queryList.begin(); p != queryList.end(); ++p)
            {
                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                positive += bf.query(*p);
                std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
                elapsedSecs += time_span.count();
            }


            std::clog << fpr << "\t" << n << "\t" << elapsedSecs << "\t" << double(positive) / queryCount - 0.5 << "\n";
        }
}



void benchmark_bf_all_positive(const uint64_t elemsStart, const uint64_t elemsEnd, const uint64_t elemsStep,
                                const double fprStart, const double fprEnd, const double fprStep,
                                const uint64_t queryCount)
{
    std::mt19937_64 rng(std::chrono::steady_clock::now().time_since_epoch().count());

    std::vector<uint64_t> keyList;


    for(double fpr = fprStart; fpr <= fprEnd; fpr += fprStep)
        for(uint64_t n = elemsStart; n <= elemsEnd; n += elemsStep)
        {
            keyList.clear();

            generate_random_list(n, keyList);

            bloom_filter bf(n, fpr);

            for(auto p = keyList.begin(); p != keyList.end(); ++p)
                bf.insert(*p);


            double elapsedSecs = 0;
            uint64_t positive = 0;
            
            for(uint64_t i = 0; i < queryCount; ++i)
            {
                uint64_t key = keyList[std::uniform_int_distribution<uint64_t>(0, n - 1)(rng)];
                std::chrono::high_resolution_clock::time_point t_start = std::chrono::high_resolution_clock::now();
                positive += bf.query(key);
                std::chrono::high_resolution_clock::time_point t_end = std::chrono::high_resolution_clock::now();

                std::chrono::duration<double> time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t_end - t_start);
                elapsedSecs += time_span.count();
            }


            std::clog << fpr << "\t" << n << "\t" << elapsedSecs << "\t" << double(positive) / queryCount - 1 << "\n";
        }
}



int main()
{
    const uint64_t elemsStart = 100000;     // 0.1M
    const uint64_t elemsEnd = 10000000;     // 10M
    const uint64_t elemsStep = 100000;      // 0.1M
    const double fprStart = 0.05;
    const double fprEnd = 0.25;
    const double fprStep = 0.05;
    const uint64_t queryCount = 100000;     // 0.1M
    
    // benchmark_bf_all_negative(elemsStart, elemsEnd, elemsStep, fprStart, fprEnd, fprStep, queryCount);

    // benchmark_bf_half_negative(elemsStep, elemsEnd, elemsStep, fprStart, fprEnd, fprStep, queryCount);

    benchmark_bf_all_positive(elemsStep, elemsEnd, elemsStep, fprStart, fprEnd, fprStep, queryCount);

    return 0;
}