#include <iostream>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <vector>

#include "murmur3.c"

class bloom_filter
{
private:
    uint64_t n;     // Expected number of distinct keys to be inserted.
    double p;       // Desired false positive rate.
    uint64_t m;     // Size of the Bloom filter.
    uint32_t k;     // Number of hash functions.

    uint64_t hashCombineFactor = 0; // A modular factor to combine the 128-bit murmur hash value
                                    // to a 64-bit one.

    std::vector<bool> B;    // The bit-array for the bloom filter.


    // Sets the hash bit-vector, along with other fields required for the Murmur3 hash.
    void set_hash_data_structure();

    // Given the expected number of distinct elements and the desired false positive
    // rate, infers the optimal setting for the bit-vector size and the number of
    // hash functions.
    void set_parameters(uint64_t numDistinctKeys, double falsePositiveRate);

    // Logs the hash bit-vector.
    void dump_bits();

    // Gets the bit index to be set / queried for a key, using the Murmur3 hash,
    // with a specific provided seed.
    inline uint64_t get_murmur_hash_slot(uint64_t &key, uint32_t &seed);


public:
    bloom_filter() {}

    // Constructs a Bloom filter with 'numDistinctKeys' expected number of distinct
    // keys, and 'falsePositiveRate' desired FPR.
    bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate);

    // Constructs a Bloom filter with 'numDistinctKeys' expected number of distinct
    // keys, and 'falsePositiveRate' desired FPR; fills it with the keys from the file
    // named 'keyFile'; and serializes the data structure to a file named 'outputFile'.
    bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile);

    // Constructs a Bloom filter from a disk-saved (serialized) file named 'bfFile';
    // i.e. deserializes the data structure from disk to memory.
    bloom_filter(std::string &bfFile);

    // Logs the parameters of the Bloom filter.
    void dump_metadata(bool dumpBits = false);

    // Inserts the key 'key' into the Bloom filter.
    void insert(uint64_t key);

    // Queries for the existence of the key 'key' into the Bloom filter.
    bool query(uint64_t key);

    // Queries for the existence of all the keys from the file named 'queryFile',
    // and put the results into the vector 'result'.
    void query(std::string &queryFile, std::vector<bool> &result);

    // Serializes the Bloom filter to the file named 'outputputFile'.
    void serialize(std::string &outputFile);

    // Deserializes a Bloom filter from the file named 'bfFile'.
    void deserialize(std::string &bfFile);

    ~bloom_filter();
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
        std::clog << B[i];

    std::clog << "\n";
}



void bloom_filter::set_parameters(uint64_t numDistinctKeys, double falsePositiveRate)
{
    n = numDistinctKeys;
    p = falsePositiveRate;

    m = ceil((-double(n) * log(p)) / (log(2) * log(2)));
    k = ceil((double(m) / n) * log(2));

    set_hash_data_structure();
}



void bloom_filter::set_hash_data_structure()
{
    // Set the bit-vector size a little larger than required, (up-to 7 bits);
    // to accomodate full byte read/write operations.
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

    input.close();

    
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

    // return ((((H[1] % m) * hashCombineFactor) % m) + (H[0] % m)) % m;
    return H[0] % m;
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



bloom_filter::~bloom_filter()
{
    if(!B.empty())
    {
        B.clear();
        B.shrink_to_fit();
    }
}
