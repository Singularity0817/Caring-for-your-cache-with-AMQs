#include <unistd.h>

#include "bloom_filter.h"

class blocked_bloom_filter
{
private:
    uint64_t n;     // Expected number of distinct keys to be inserted.
    double p;       // Desired false positive rate.
    uint64_t m;     // Size of the bloom filter.
    uint32_t k;     // Number of hash functions.

    uint64_t hashCombineFactor = 0; // A modular factor to combine the 128-bit murmur hash value
                                    // to a 64-bit one.
    
    uint64_t blockCount;    // Number of blocks (bloom filters).
    uint64_t blockSize;     // Size of each block (bloom filter), in bits.

    std::vector<bloom_filter> BF;   // Vector of the Bloom filter blocks.


    // Set the expected number of distinct elements and the desired false positive rate;
    // infer the optimal setting for the bit-vector size and the number of hash functions.
    // Also, set the block size for each Bloom filter block to the cache line size.
    void set_parameters(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSizeinBits = 0);

    // Get the block number of the Bloom filter to be modified / queried for a key,
    // using the Murmur3 hash.
    inline uint64_t get_murmur_hash_block(uint64_t &key);

    // Serialize the blocked Bloom filter to the file named 'outputputFile'.
    void serialize(std::string &outputFile);

    // Deserialize a Bloom filter from the file named 'bfFile'.
    void deserialize(std::string &bbfFile);


public:

    blocked_bloom_filter() {}

    // Construct a blocked Bloom filter with 'numDistinctKeys' expected number of
    // distinct keys, and 'falsePositiveRate' desired FPR. Also, set the block size
    // for each Bloom filter block to the cache line size; if you are not familiar
    // with the cache line staff, you can (and are encouraged to) leave using the
    // last parameter.
    blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSizeinBits = 0);

    // Construct a blocked Bloom filter with numDistinctKeys' expected number of
    // distinct keys, and 'falsePositiveRate' desired FPR; fill it with the keys
    // from the file named 'keyFile', and serialize the data structure to a file
    // named 'outputFile'.
    blocked_bloom_filter(std::string &keyFile, uint64_t numDistinct, double fpr, std::string &outputFile);

    // Construct a Bloom filter from a disk-saved (serialized) file named 'bbfFile';
    // i.e. deserialize the data structure from disk to memory.
    blocked_bloom_filter(std::string &bbfFile);

    // Log the parameters of the blocked Bloom filter.
    void dump_metadata();

    // Insert the key 'key' into the blocked Bloom filter.
    void insert(uint64_t key);

    // Query for the existence of the key 'key' into the blocked Bloom filter.
    bool query(uint64_t key);
};



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



void blocked_bloom_filter::set_parameters(uint64_t numDistinctKeys, double falsePositiveRate, uint64_t blockSizeinBits)
{
    n = numDistinctKeys;
    p = falsePositiveRate;
    blockSize = (blockSizeinBits ? blockSizeinBits : 8 * sysconf(_SC_LEVEL1_DCACHE_LINESIZE));

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
        

    hashCombineFactor = (((uint64_t(1) << 63) % blockCount) * (2 % blockCount)) % blockCount;
}



blocked_bloom_filter::blocked_bloom_filter(uint64_t numDistinctKeys, double falsePositiveRate,
                                            uint64_t blockSizeinBits)
{
    set_parameters(numDistinctKeys, falsePositiveRate, blockSizeinBits);

    dump_metadata();
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

    MurmurHash3_x64_128(&key, sizeof(key), k - 1, H);

    return ((((H[1] % blockCount) * hashCombineFactor) % blockCount) + (H[0] % blockCount)) % blockCount;
}



void blocked_bloom_filter::insert(uint64_t key)
{
    uint64_t block = get_murmur_hash_block(key);

    BF[block].insert(key);
}



bool blocked_bloom_filter::query(uint64_t key)
{
    uint64_t block = get_murmur_hash_block(key);

    return BF[block].query(key);
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


    // Serialize the blocks (the Bloom filters).

    for(uint64_t i = 0; i < blockCount; ++i)
        BF[i].serialize(output, true);

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

    set_parameters(n, p, blockSize);


    // Deserialize the blocks (the Bloom filters).

    for(uint64_t i = 0; i < blockCount; ++i)
        BF[i].deserialize(input, true);

    input.close();
}
