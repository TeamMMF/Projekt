//
// Created by matep on 08/12/2017.
//

typedef struct {
    uint64_t hash;
    uint16_t index;
    bool rev;
} minimizer;

typedef struct {
    uint16_t index;
    bool rev;
} hashEntry;

typedef struct {
    uint64_t hash;
    char* seq_id;
    uint16_t index;
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
