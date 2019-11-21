#include <unistd.h>

#include "bloom_filter.h"

class blocked_bloom_filter
{
private:
    uint64_t n;     // Expected number of distinct keys to be inserted.
    double p;       // Desired false positive rate.
    uint64_t m;     // Size of the bloom filter.
    uint32_t k;     // Number of hash functions.
    
    uint64_t blockCount;    // Number of blocks (bloom filters).
    uint64_t blockSize;     // Size of each block (bloom filter), in bits.

    std::vector<bloom_filter> BF;


    void set_parameters(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSize);


public:

    blocked_bloom_filter() {}

    blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSize = 0);

    void dump_metadata();
};



blocked_bloom_filter::blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSize)
{
    set_parameters(numDistinctKeys, falsePositiveRate,
                    8 * (blockSize ? blockSize : sysconf(_SC_LEVEL1_DCACHE_LINESIZE)));

    dump_metadata();
}



void blocked_bloom_filter::set_parameters(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSize)
{
    n = numDistinctKeys;
    p = falsePositiveRate;
    this -> blockSize = blockSize;

    m = ceil((-double(n) * log(p)) / (log(2) * log(2)));
    if(m % blockSize)
        m += (blockSize - m % blockSize);

    k = (double(m) / n) * log(2);
    if(k == 1)
        k = 2;

    blockCount = m / blockSize;
    BF.reserve(blockCount);

    for(uint64_t i = 0; i < blockCount; ++i)
        BF.emplace_back(blockSize, k - 1);
}



void blocked_bloom_filter::dump_metadata()
{
    std::clog << "\tn = " << n << ";";
    std::clog << "\tp = " << p << ";";
    std::clog << "\tm = " << m << ";";
    std::clog << "\tk = " << k << ";";
    std::clog << "\n";

    std::clog << "\tFilter count = " << blockCount << ";";
    std::clog << "\tFilter size = " << blockSize << ";";
    std::clog << "\n";
}
