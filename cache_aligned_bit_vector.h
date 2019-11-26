#include <cstdint>
#include <cstdio>
#include <cstring>
#include <cstdlib>
#include <fstream>
#include <random>
#include <functional>
#include <unistd.h>


class bit_vector
{
private:
    const static int UNIT_WIDTH = 8;
    const static uint64_t ALL_ONES_UNIT = 0xFF;

    uint64_t len;
    unsigned char *B;

public:
    bit_vector() { len = 0; }
    bit_vector(uint64_t len);
    bit_vector(bool *bits, uint64_t len);
    
    ~bit_vector() { free(B); }


    inline uint64_t get_len() { return len; }
    inline void set_len(uint64_t len);
    inline bool get_bit(uint64_t idx);
    inline void set_bit(uint64_t idx);
    inline void reset_bit(uint64_t idx);
    uint64_t get_int(uint64_t idx, uint64_t len);
    void set_int(uint64_t idx, uint64_t len, uint64_t val);
    void print();
    void serialize(std::ofstream &output);
    void deserialize(std::ifstream &input);
    void generate_random_bitvector(uint64_t length);
};



bit_vector::bit_vector(uint64_t len)
{
    this -> len = len;

    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;
    B = (unsigned char *) aligned_alloc(sysconf(_SC_LEVEL1_DCACHE_LINESIZE), slots * sizeof(unsigned char));
}



bit_vector::bit_vector(bool *bits, uint64_t len)
{
    this -> len = len;
    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;
    B = (unsigned char *) aligned_alloc(sysconf(_SC_LEVEL1_DCACHE_LINESIZE), slots * sizeof(unsigned char));

    for(uint64_t i = 0; i < len; ++i)
        bits[i] ? set_bit(i) : reset_bit(i);
}



void bit_vector::set_len(uint64_t len)
{
    if(this -> len)
        free(B);

    this -> len = len;

    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;
    B = (unsigned char *) aligned_alloc(sysconf(_SC_LEVEL1_DCACHE_LINESIZE), slots * sizeof(unsigned char));

    memset(B, 0, slots * sizeof(unsigned char));
}



bool bit_vector::get_bit(uint64_t idx)
{
    return (B[idx / UNIT_WIDTH]) & (1u << (idx % UNIT_WIDTH));
}



void bit_vector::set_bit(uint64_t idx)
{
    B[idx / UNIT_WIDTH] |= (1u << (idx % UNIT_WIDTH));
}



void bit_vector::reset_bit(uint64_t idx)
{
    B[idx / UNIT_WIDTH] &= ~(1u << (idx % UNIT_WIDTH));
}



void bit_vector::set_int(uint64_t idx, uint64_t len, uint64_t val)
{
    while(len && idx % UNIT_WIDTH)
    {
        (val & 1) ? set_bit(idx) : reset_bit(idx);
        idx++, len--, val >>= 1;
    }

    while(len >= UNIT_WIDTH)
    {
        B[idx / UNIT_WIDTH] = (val & ALL_ONES_UNIT);
        idx += UNIT_WIDTH, len -= UNIT_WIDTH, val >>= UNIT_WIDTH;
    }

    while(len)
    {
        (val & 1) ? set_bit(idx) : reset_bit(idx);
        idx++, len--, val >>= 1;
    }
}



uint64_t bit_vector::get_int(uint64_t idx, uint64_t len)
{
    uint64_t val = 0, i = 0;

    while(i < len && idx % UNIT_WIDTH)
    {
        if(get_bit(idx))
            val |= (1 << i);

        idx++, i++;
    }

    while((len - i) >= UNIT_WIDTH)
    {
        val |= (B[idx / UNIT_WIDTH] << i);
        idx += UNIT_WIDTH, i += UNIT_WIDTH;
    }

    while(i < len)
    {
        if(get_bit(idx))
            val |= (1 << i);

        idx++, i++;
    }

    return val;
}



void bit_vector::print()
{
    for(uint16_t i = 0; i < len; ++i)
        putchar(get_bit(i) ? '1' : '0');
    putchar('\n');
}



void bit_vector::serialize(std::ofstream &output)
{
    output.write((const char *)&len, sizeof(len));
    
    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;

    output.write((const char *)B, slots);
}



void bit_vector::deserialize(std::ifstream &input)
{
    uint64_t l;
    input.read((char *)&l, sizeof(l));
    
    set_len(l);
    
    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;
    input.read((char *)B, slots);
}



void bit_vector::generate_random_bitvector(uint64_t length)
{
    if(this -> len)
        free(B);

    this -> len = length;
    uint64_t slots = (len + UNIT_WIDTH - 1) / UNIT_WIDTH;
    B = (unsigned char *) aligned_alloc(sysconf(_SC_LEVEL1_DCACHE_LINESIZE), slots * sizeof(unsigned char));

    
    auto gen = std::bind(std::uniform_int_distribution<>(0,1), std::default_random_engine());
    for(uint64_t idx = 0; idx < len; ++idx)
        gen() ? set_bit(idx) : reset_bit(idx);
}