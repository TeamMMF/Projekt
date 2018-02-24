//
// Created by matep on 08/12/2017.
//

#ifndef CustomTypes
#define CustomTypes

#include <cstdint>
/**
 * Represents a minimizer tuple (h,i) where
 * h - hash
 * i - start index in the sequence (positive if minimizer is not from reversed starnd, negative otherwise)
 */
typedef struct {
    uint64_t hash;
    int32_t index;
} minim;

/**
 * Represents a tuple which carries information about the sequece and start index in that seqeunce for a certain minimizer
 * hash value.
 */
typedef struct {
    uint32_t seq_id;
    int32_t index;
} hashMinPair3;

#endif
