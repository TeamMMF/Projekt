#include "ProjectConfig.h"

#include <string>
#include <sstream>
#include <iostream>
#include <vector>
#include <memory>
#include <stdint.h>
#include <string.h>
#include "bioparser/bioparser.hpp"

using namespace std;

// define a class for sequences in FASTA format
class FASTASampleClas {

	char * name;
	char * data;
	string test_string ="asdfasdfasfd";

	uint32_t name_length;
	uint32_t data_length;

    public:
        // required signature for the constructor
        FASTASampleClas(uint64_t object_id,
            const char* name,
            uint32_t name_length,
            const char* data,
            uint32_t data_length) {

        	   this->name_length = name_length;
        	   this->name =  new char[name_length];
        	   strncpy(this->name,name,name_length);

        	   this->data_length = data_length;
        	   this->data=new char[data_length];
        	   strncpy(this->data,data,data_length);
        }

        string get_description();
        string get_test_string();
};

string FASTASampleClas :: get_description(){
	ostringstream oss;
	
	oss << "Name: ";
	for (int i = 0; i < name_length; ++i)
	{
		oss << name[i];
	}

	oss << endl;
	oss << "Data: ";
	for (int i = 0; i < data_length; ++i)
	{
		oss << data[i];
	}

	return oss.str();
}

string FASTASampleClas :: get_test_string(){
	return test_string;
}

int main(int argc, char const *argv[])
{
	vector<unique_ptr<FASTASampleClas>> fasta_objects;

	string project_root(PROJECT_ROOT);
	auto fasta_reader = bioparser::createReader<FASTASampleClas, bioparser::FastaReader>(project_root+"src/resources/sample.fasta");
	
	fasta_reader->read_objects(fasta_objects, -1);

	for (vector<unique_ptr<FASTASampleClas>>::iterator i = fasta_objects.begin(); i != fasta_objects.end(); ++i)
	{
		cout << (**i).get_description() << endl << endl;
	}

	return 0;
}
