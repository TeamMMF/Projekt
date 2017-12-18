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
    int target_min_index;
    int target_max_index;
    int query_min_index;
    int query_max_index;
    //string query_name;
    bool reverse;
} mapInfo;

typedef struct {
    uint32_t index;
    bool rev;
} hashEntry;


typedef struct {
    uint64_t hash;
    uint32_t seq_id;
    uint32_t index;
    bool rev;
} hashMinPair;

typedef struct {
    uint64_t hash;
    const char* seq_id;
    uint32_t seq_id_length;
    uint32_t index;
    bool rev;
} hashMinPair2;

typedef struct{
    char* qs_name;
    int qs_length;
    int q_start;
    int q_end;
    char strand;
    char* ts_name;
    int ts_length;
    int target_start;
    int target_end;
    int residue;
    int alignment_block_end;
    int mapping_quality;
} PAF_data;
#endif
