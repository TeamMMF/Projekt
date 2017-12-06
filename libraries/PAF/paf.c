//
// Created by filip on 05.12.17..
//


typedef PAF_data struct{
    char* qs_name;
    int qs_length;
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
};

public int write_paf_to_file(char *output_file, PAF_data *data, int size){
    FILE *file = fopen(output_file, "w+");
    for (int i = 0; i < size; ++i) {
        fprintf(file, "%s\t%d\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d\n",
                data[i].qs_name,
                data[i].qs_length,
                data[i].q_start,
                data[i].q_end,
                data[i].strand,
                data[i].ts_name,
                data[i].ts_length,
                data[i].target_start,
                data[i].target_end,
                data[i].residue,
                data[i].alignment_block_end,
                data[i].mapping_quality
        );
    }
    fclose(file);
}

