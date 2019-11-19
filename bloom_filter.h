#include <cstdint>
#include <cmath>
#include <vector>
#include <iostream>

#include "murmur3.c"

class bloom_filter
{
private:
    uint64_t n;     // Expected number of distinct keys to be inserted.
    double p;       // Desired false positive rate.
    uint64_t m;     // Size of the bloom filter.
    uint32_t k;     // Number of hash functions.

    uint64_t hashCombineFactor = 0; // A modular factor to combine the 128-bit murmur hash value
                                    // to a 64-bit one.

    std::vector<bool> B;    // The bit-array for the bloom filter.


    inline uint64_t get_murmur_hash_slot(uint64_t &key, uint32_t &seed);

public:
    bloom_filter() {}

    bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate);

    void dump_metadata();

    void insert(uint64_t key);

    bool query(uint64_t key);
};



bloom_filter::bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate):
    n (numDistinctKeys),
    p (falsePositiveRate)
{
    B.resize(n);

    m = ceil((-double(n) * log(p)) / (log(2) * log(2)));
    k = (double(m) / n) * log(2);

    hashCombineFactor = (((uint64_t(1) << 63) % m) * (2 % m)) % m;
}



void bloom_filter::dump_metadata()
{
    std::clog << "n = " << n << "; p = " << p << "; m = " << m << "; hc = " << hashCombineFactor << "\n";
}



uint64_t bloom_filter::get_murmur_hash_slot(uint64_t &key, uint32_t &seed)
{
    uint64_t H[2];

    MurmurHash3_x64_128(&key, sizeof(key), seed, H);

    return ((((H[1] % m) * hashCombineFactor) % m) + (H[0] % m)) % m;
}



void bloom_filter::insert(uint64_t key)
{
    for(uint32_t seed = 0; seed < k; ++seed)
        B[get_murmur_hash_slot(key, seed)] = true;
}



bool bloom_filter::query(uint64_t key)
{
    for(uint32_t seed = 0; seed < k; ++seed)
        if(!B[get_murmur_hash_slot(key, seed)])
            return false;

    return true;
}