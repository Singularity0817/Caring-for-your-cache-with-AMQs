# Caring for Your Cache with AMQs


Overview
--------

Implementation of: an efficient basic Bloom filter library and a variant of it, the blocked Bloom filter (with cache-alignment support) library; with insert, query, and (de)serialization support.

Compile
--------
```
g++ -std=c++17 -o bf bf.cpp
g++ -std=c++17 -o bbf bbf.cpp
```

API
--------
* `./bf build -k <key file> -f <fpr> -n <num. distinct keys> -o <output file>`: reads in the keys in the file `<key file>`, and constructs a standard Bloom filter with a target false positive rate of `<fpr>`, and then serializes the constructed Bloom filter to the file `<output file>`. The `<key file>` contains a linebreak-separated list of keys. The `<num. distinct keys>` parameter is required to set the internal hashing data structure parameters of the Bloom filter; and should be set to a good estimate of the number of distinct keys present at the `<key file>`, or higher than that if memory-consumption is not much of an issue.

* `./bf query -i <input file> -q <queries>`: Loads a serialized standard Bloom filter from the file `<input file>`, and issues a series of queries present in the file `<queries>`. The `<queries>` file contains linebreak-separated list of keys to be queried. The query results are displayed to the standard out.

* `./bbf build -k <key file> -f <fpr> -n <num. distinct keys> -o <output file>`: reads in the keys in the file `<key file>`, and constructs a blocked Bloom filter with a target false positive rate of `<fpr>`, and then serializes the constructed blocked Bloom filter to the file `<output file>`. The `<key file>` contains a linebreak-separated list of keys. The `<num. distinct keys>` parameter is required to set the internal hashing data structure parameters of the blocked Bloom filter; and should be set to a good estimate of the number of distinct keys present at the `<key file>`, or higher than that if memory-consumption is not much of an issue.

* `./bbf query -i <input file> -q <queries>`: Loads a serialized blocked Bloom filter from the file `<input file>`, and issues a series of queries present in the file `<queries>`. The `<queries>` file contains linebreak-separated list of keys to be queried. The query results are displayed to the standard out.
