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



    void set_parameters(uint64_t numDistinctKeys, double falsePositiveRate);

    inline uint64_t get_murmur_hash_slot(uint64_t &key, uint32_t &seed);

    void serialize(std::string &outputFile);

public:
    bloom_filter() {}

    bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate);

    bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile);

    void dump_metadata();

    void insert(uint64_t key);

    bool query(uint64_t key);
};



void bloom_filter::dump_metadata()
{
    std::clog << "\tn = " << n << ";";
    std::clog << "\tp = " << p << ";";
    std::clog << "\tm = " << m << ";";
    std::clog << "\tk = " << k << ";";
    std::clog << "\thc = " << hashCombineFactor << ";";
    std::clog << "\n";
}



void bloom_filter::set_parameters(uint64_t numDistinctKeys, double falsePositiveRate)
{
    n = numDistinctKeys;
    p = falsePositiveRate;

    m = ceil((-double(n) * log(p)) / (log(2) * log(2)));
    k = (double(m) / n) * log(2);

    B.resize(m + (m % 8 ? (8 - m % 8) : 0));

    hashCombineFactor = (((uint64_t(1) << 63) % m) * (2 % m)) % m;
}



bloom_filter::bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate)
{
    set_parameters(numDistinctKeys, falsePositiveRate);
}



bloom_filter::bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile)
{
    std::ifstream input(keyFile);
    if(!input.is_open())
    {
        std::cerr << "Cannot open keys file " << keyFile << "\n";
        exit(1);
    }

    set_parameters(numDistinct, fpr);
    dump_metadata();

    uint64_t key;
    while(input >> key)
        insert(key);

    
    serialize(outputFile);
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



void bloom_filter::serialize(std::string &outputFile)
{
    std::ofstream output;
    output.open(outputFile.c_str(), std::ios::binary | std::ios::out);


    // Serialize the parameters.
    
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)&p, sizeof(p));
    output.write((const char *)&m, sizeof(m));
    output.write((const char *)&k, sizeof(k));


    // Serialize the bits array.

    for(uint64_t i = 0; i < m; i += 8)
    {
        unsigned char byte = 0;
        for(int j = 0; j < 8; ++j)
            byte |= B[i + j];

        output.write((const char *)&byte, sizeof(byte));
    }


    output.close();
}
