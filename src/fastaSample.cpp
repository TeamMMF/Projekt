#include "ProjectConfig.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>

#include "FASTASampleClass.cpp"
#include "bioparser/bioparser.hpp"

using namespace std;

int main(int argc, char const *argv[])
{

	string project_root(PROJECT_ROOT);

    vector<unique_ptr<FASTASampleClass>> fasta_objects;
	auto fasta_reader = bioparser::createReader<FASTASampleClass, bioparser::FastaReader>(project_root + "src/resources/sample.fasta");
	
	fasta_reader->read_objects(fasta_objects, static_cast<uint64_t>(-1));

	for (auto &fasta_object : fasta_objects) {
		cout << (*fasta_object).get_description() << endl << endl;
	}

	return 0;
}
