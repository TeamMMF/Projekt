#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <stdlib.h>

#define MAX_OCTETS 256

typedef unsigned char ubyte;

struct t{
    short is_character;
    ubyte value;
    struct t* left;
    struct t* right;
};
typedef struct t node;

FILE *safe_fopen(char* filepath, char* option){
    FILE *file = fopen(filepath, option);
    if(file == NULL){
        fprintf(stderr,"Could not open file '%s' with option '%s'.\n", filepath, option);
        exit(1);
    }
}

node *new_node(short is_character, ubyte value) {
    node *new_node = (node *) malloc(sizeof(node));

    new_node->is_character = is_character;
    new_node->value = value;
    new_node->right = NULL;
    new_node->left = NULL;

    return new_node;
}

void read_table_file(char *input_file, ubyte data[MAX_OCTETS][MAX_OCTETS]) {
    FILE *file = safe_fopen(input_file, "r");
    int i = 0;
    while(i<MAX_OCTETS){
        if(!fgets(data[i],MAX_OCTETS,file)){
            fprintf(stderr,"Invalid table file");
            exit(1);
        }
        strtok(data[i++],"\n");
    }
    fclose(file);
}

/*
void add_to_tree(ubyte value, ubyte* code, int index, node *current){
    if(current == NULL){
        current = new_node(0,0);
    }
    if(!code[index]){
        current -> is_character = 1;
        current -> value = value;
        return;
    }
    printf("%d\n",value);
    if(code[index] == '1'){
        root->right = add_to_tree(value, code, index + 1, next);
    }else{
        root->left = add_to_tree(value, code, index + 1, current -> left);
    }
    node* next = code[index] == '1' ? current -> right : current -> left;
    add_to_tree(value, code, index + 1, next);
}
*/

void add_to_tree(ubyte value, ubyte code[],node *root){
    int index = 0;
    while(code[index]!='\0'){
        switch(code[index]){
            case '0':
                if(root -> left == NULL){
                    root -> left = new_node(0,0);
                }
                root = root -> left;
                break;
            case '1':
                if(root -> right == NULL){
                    root -> right = new_node(0,0);
                }
                root = root -> right;
                break;
            default:
                fprintf(stderr, "Illegal symbol '%c' found in the table file.\n", code[index]);
                exit(1);
        }
        index++;
    }
    root -> is_character = 1;
    root -> value = value;
}

node *construct_decision_tree(ubyte data[MAX_OCTETS][MAX_OCTETS]){
    node *root = new_node(0,0);
    for (int i = 0; i < MAX_OCTETS; ++i)
    {
        add_to_tree(i, data[i], root);
    }
    return root;
}

node* transition(node* start, short symbol){
    if(symbol == 1){
        return start -> right;
    }else if( symbol == 0){
        return start -> left;
    }else{
        return NULL;
    }
}

node *decode_file(char* input_file, char* output_file, node* decision_root){
    FILE *in = safe_fopen(input_file, "rb");
    FILE *out = safe_fopen(output_file, "w");

    ubyte in_buff = 0;
    int in_buff_pos = 7;

    node* current = decision_root;
    while(fread(&in_buff,1,1,in) == 1){
        short is_char;
        do{
            current = transition(current,(in_buff>>in_buff_pos) & 1);
            is_char = current -> is_character;
            in_buff_pos--;

            if(is_char){
                fprintf(out,"%c", current -> value);
                current = decision_root;
            }

        }while(in_buff_pos != -1);

        in_buff = 0;
        in_buff_pos = 7;
    }

    fclose(in);
    fclose(out);
}

void destroy(node* root){
    if(!root){
        return;
    }

    destroy(root->left);
    destroy(root->right);

    free(root);
}

int main(int argc, char *argv[])
{
    ubyte data[MAX_OCTETS][MAX_OCTETS];
    read_table_file(argv[1], data);
    node* decision_tree = construct_decision_tree(data);
    decode_file(argv[2], argv[3], decision_tree);
    destroy(decision_tree);
    printf("Successfully decoded file '%s' into file '%s'.\n", argv[2], argv[3]);
    return 0;
}