#include <iostream>
#include <fstream>
#include <cstdint>
#include <cmath>
#include <vector>
#include <unistd.h>

#include "murmur3.c"


class blocked_bloom_filter
{
private:
    uint64_t n;     // Expected number of distinct keys to be inserted.
    double p;       // Desired false positive rate.
    uint64_t m;     // Size of the blocked Bloom filter.
    uint32_t k;     // Number of hash functions.

    // Modular factors to combine the 128-bit murmur hash values to 64-bit ones.
    uint64_t hashCombineFactor_block; 
    uint64_t hashCombineFactor_slot;
    
    uint64_t blockCount;    // Number of filter blocks.
    uint64_t blockSize;     // Size of each filter block in bits.

    std::vector<bool> B;    // The bit-vector for the blocked Bloom filter.


    // Sets the hash bit-vector, along with other fields required for the
    // Murmur3 hash.
    void set_hash_data_structure();

    // Given the expected number of distinct elements and the desired false positive
    // rate, infers the optimal setting for the bit-vector size and the number of hash
    // functions. Also, sets the block size for each block to the cache line size
    // (or the provided size).
    void set_parameters(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSizeinByte = 0);

    // Get the block number of the filter block to be modified / queried for a key,
    // using the Murmur3 hash.
    inline uint64_t get_murmur_hash_block(uint64_t &key);

    // Get the bit index to be set / queried for a key, using the Murmur3 hash,
    // with a specific provided seed.
    inline uint64_t get_murmur_hash_slot(uint64_t &key, uint32_t &seed);

    // Serialize the blocked Bloom filter to the file named 'outputputFile'.
    void serialize(std::string &outputFile);

    // Deserialize a Bloom filter from the file named 'bbfFile'.
    void deserialize(std::string &bbfFile);


public:

    blocked_bloom_filter() {}

    // Constructs a blocked Bloom filter with 'numDistinctKeys' expected number of
    // distinct keys, and 'falsePositiveRate' desired FPR. Also, sets the block size
    // for each filter block to the cache line size; if you are not familiar
    // with the cache line staff, you can (and are encouraged to) leave using the
    // last parameter.
    blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSizeinByte = 0);

    // Constructs a blocked Bloom filter with 'numDistinctKeys' expected number of
    // distinct keys, and 'falsePositiveRate' desired FPR; fills it up with the keys
    // from the file named 'keyFile'; and serializes the data structure to a file
    // named 'outputFile'.
    blocked_bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile);

    // Constructs a blocked Bloom filter from a disk-saved (serialized) file named
    // 'bbfFile'; i.e. deserializes the data structure from disk to memory.
    blocked_bloom_filter(std::string &bbfFile);

    // Logs the parameters of the blocked Bloom filter.
    void dump_metadata();

    // Inserts the key 'key' into the blocked Bloom filter.
    void insert(uint64_t key);

    // Queries for the existence of the key 'key' into the blocked Bloom filter.
    bool query(uint64_t key);

    // Queries for the existence of all the keys from the file named 'queryFile',
    // and puts the results into the vector 'result'.
    void query(std::string &queryFile, std::vector<bool> &result);

    ~blocked_bloom_filter();
};



void blocked_bloom_filter::dump_metadata()
{
    std::clog << "\tn = " << n << ";";
    std::clog << "\tp = " << p << ";";
    std::clog << "\tm = " << m << ";";
    std::clog << "\tk = " << k << ";";
    std::clog << "\n";

    std::clog << "\tBlock count = " << blockCount << ";";
    std::clog << "\tBlock size (in bits) = " << blockSize << ";";
    std::clog << "\n";
}



void blocked_bloom_filter::set_parameters(uint64_t numDistinctKeys, double falsePositiveRate,
                                            uint64_t blockSizeinByte)
{
    n = numDistinctKeys;
    p = falsePositiveRate;
    blockSize = 8 * (blockSizeinByte ? blockSizeinByte : sysconf(_SC_LEVEL1_DCACHE_LINESIZE));

    m = ceil((-double(n) * log(p)) / (log(2) * log(2)));
    if(m % blockSize)
        m += (blockSize - m % blockSize);

    k = ceil((double(m) / n) * log(2));
    if(k == 1)
        k = 2;

    blockCount = m / blockSize;


    set_hash_data_structure();
}



void blocked_bloom_filter::set_hash_data_structure()
{
    // Sets the bit-vector size a little larger than required, (up-to 7 bits);
    // to accomodate full byte read/write operations.
    B.resize(m + (m % 8 ? (8 - m % 8) : 0));

    hashCombineFactor_block = (((uint64_t(1) << 63) % blockCount) * (2 % blockCount)) % blockCount;
    hashCombineFactor_slot = (((uint64_t(1) << 63) % blockSize) * (2 % blockSize)) % blockSize;
}



blocked_bloom_filter::blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate,
                                            uint64_t blockSizeinByte)
{
    set_parameters(numDistinctKeys, falsePositiveRate, blockSizeinByte);

    // dump_metadata();
}


blocked_bloom_filter::blocked_bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr,
                                            std::string &outputFile)
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



blocked_bloom_filter::blocked_bloom_filter(std::string &bbfFile)
{
    deserialize(bbfFile);

    dump_metadata();
}



uint64_t blocked_bloom_filter::get_murmur_hash_block(uint64_t &key)
{
    uint64_t H[2];

    MurmurHash3_x64_128(&key, sizeof(key), 0, H);

    // return ((((H[1] % blockCount) * hashCombineFactor_block) % blockCount) + (H[0] % blockCount)) % blockCount;
    return H[0] % blockCount;
}



uint64_t blocked_bloom_filter::get_murmur_hash_slot(uint64_t &key, uint32_t &seed)
{
    uint64_t H[2];

    MurmurHash3_x64_128(&key, sizeof(key), seed, H);

    // return ((((H[1] % blockSize) * hashCombineFactor_slot) % blockSize) + (H[0] % blockSize)) % blockSize;
    return H[0] % blockCount;
}



void blocked_bloom_filter::insert(uint64_t key)
{
    uint64_t block = get_murmur_hash_block(key);

    uint64_t offset = block * blockSize;
    for(uint32_t seed = 1; seed < k; ++seed)
        B[offset + get_murmur_hash_slot(key, seed)] = true;
}



bool blocked_bloom_filter::query(uint64_t key)
{
    uint64_t block = get_murmur_hash_block(key);

    uint64_t offset = block * blockSize;
    for(uint32_t seed = 1; seed < k; ++seed)
        if(!B[offset + get_murmur_hash_slot(key, seed)])
            return false;

    return true;
}



void blocked_bloom_filter::query(std::string &queryFile, std::vector<bool> &result)
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



void blocked_bloom_filter::serialize(std::string &outputFile)
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
    output.write((const char *)&blockSize, sizeof(blockSize));


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



void blocked_bloom_filter::deserialize(std::string &bbfFile)
{
    std::ifstream input;

    input.open(bbfFile.c_str(), std::ios::binary | std::ios::in);

    if(!input.is_open())
    {
        std::cerr << "Cannot open output file " << bbfFile << "\n";
        exit(1);
    }


    // Deserialize the parameters.

    input.read((char *)&n, sizeof(n));
    input.read((char *)&p, sizeof(p));
    input.read((char *)&blockSize, sizeof(blockSize));

    set_parameters(n, p, blockSize / 8);


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



blocked_bloom_filter::~blocked_bloom_filter()
{
    if(!B.empty())
    {
        B.clear();
        B.shrink_to_fit();
    }
}
