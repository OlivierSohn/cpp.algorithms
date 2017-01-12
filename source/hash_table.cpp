
namespace imajuscule {

    // good_primes is made according to
    // http://planetmath.org/goodhashtableprimes :
    //
    //   - 2^(i+1) >= good_primes_for_hash_function[i] > 2^i
    //   - the distance between prime and nearest powers of 2 is maximal
    
    UniversalHashFunction::primes UniversalHashFunction::good_primes {
        {
            2,
            3,
            5,
            13,
            23,
            53,
            97,
            193,
            389,
            769,
            1543,
            3079,
            6151,
            12289,
            24593,
            49157,
            98317,
            196613,
            393241,
            786433,
            1572869,
            3145739,
            6291469,
            12582917,
            25165843,
            50331653,
            100663319,
            201326611,
            402653189,
            805306457,
            1610612741
        }
    };
    
} // NS imajuscule

