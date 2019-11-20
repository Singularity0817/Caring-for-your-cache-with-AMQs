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

    void dump_bits();

    inline uint64_t get_murmur_hash_slot(uint64_t &key, uint32_t &seed);

    void serialize(std::string &outputFile);

    void deserialize(std::string &bfFile);

public:
    bloom_filter() {}

    bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate);

    bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile);

    bloom_filter(std::string &bfFile);

    void dump_metadata(bool dumpBits = false);

    void insert(uint64_t key);

    bool query(uint64_t key);

    void query(std::string &queryFile, std::vector<bool> &result);
};



void bloom_filter::dump_metadata(bool dumpBits)
{
    std::clog << "\tn = " << n << ";";
    std::clog << "\tp = " << p << ";";
    std::clog << "\tm = " << m << ";";
    std::clog << "\tk = " << k << ";";
    std::clog << "\thc = " << hashCombineFactor << ";";
    std::clog << "\n";


    if(dumpBits)
        dump_bits();
}



void bloom_filter::dump_bits()
{
    for(uint64_t i = 0; i < m; ++i)
        std::cout << B[i];

    std::cout << "\n";
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

    uint64_t key;
    while(input >> key)
        insert(key);

    
    serialize(outputFile);

    dump_metadata();
}



bloom_filter::bloom_filter(std::string &bfFile)
{
    deserialize(bfFile);

    dump_metadata();
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



void bloom_filter::query(std::string &queryFile, std::vector<bool> &result)
{
    std::ifstream input(queryFile);

    if(!input.is_open())
    {
        std::cerr << "Cannot open queries file " << queryFile << "\n";
        exit(1);
    }


    uint64_t key;

    uint64_t queryCount = 0;
    while(input >> key)
        result.push_back(query(key));

    input.close();
}



void bloom_filter::serialize(std::string &outputFile)
{
    std::ofstream output;
    output.open(outputFile.c_str(), std::ios::binary | std::ios::out);

    if(!output.is_open())
    {
        std::cerr << "Cannot open output file " << outputFile << "\n";
        exit(1);
    }


    // Serialize the parameters.
    
    output.write((const char *)&n, sizeof(n));
    output.write((const char *)&p, sizeof(p));


    // Serialize the bits array.

    for(uint64_t i = 0; i < m; i += 8)
    {
        unsigned char byte = 0;
        for(int j = 0; j < 8; ++j)
            byte |= (B[i + j] << j);

        output.write((const char *)&byte, sizeof(byte));
    }


    output.close();
}



void bloom_filter::deserialize(std::string &bfFile)
{
    std::ifstream input;

    input.open(bfFile.c_str(), std::ios::binary | std::ios::in);

    if(!input.is_open())
    {
        std::cerr << "Cannot open output file " << bfFile << "\n";
        exit(1);
    }


    // Deserialize the parameters.

    input.read((char *)&n, sizeof(n));
    input.read((char *)&p, sizeof(p));

    set_parameters(n, p);


    // Deserialize the bits array.

    for(uint64_t i = 0; i < m; i += 8)
    {
        unsigned char byte;
        input.read((char *)&byte, sizeof(byte));

        for(int j = 0; j < 8; ++j)
            B[i + j] = (byte & (1 << j));
    }


    input.close();
}
