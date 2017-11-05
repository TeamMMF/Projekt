#include <iostream>
#include <stdio.h>



int main(int argc, char const *argv[])
{
	auto input_file = fopen("/home/filip/Desktop/Projekt/Projekt/src/resources/sample.fasta", "r");
    if (input_file == nullptr) {
        fprintf(stderr, "bioparser::createReader error: "
            "unable to open file %s!\n", "asdf");
    }
    std::cout << errno << std::endl;

}