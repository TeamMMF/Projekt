//
// Created by matep on 08/12/2017.
//

#ifndef CustomTypes
#define CustomTypes

#include <cstdint>

typedef struct {
    uint64_t hash;
    uint32_t index;
    bool rev;
} minimizer;

typedef struct {
    uint64_t hash;
    int32_t index;
} minim;

typedef struct {
    uint32_t seq_id;
    uint32_t index;
    bool rev;
} hashMinPair2;

typedef struct {
    uint32_t seq_id;
    int32_t index;
} hashMinPair3;

#endif
