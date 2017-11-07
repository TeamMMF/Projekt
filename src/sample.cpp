#include "ProjectConfig.h"

#include <cstdio>
#include "Common.hpp"


int main(int argc, char const *argv[])
{
	if (argc < 2){
	    fprintf(stdout,"%s Version %d.%d\n",
	            argv[0],
	            SequenceOverlaping_VERSION_MAJOR,
	            SequenceOverlaping_VERSION_MINOR);
	    fprintf(stdout,"Usage: %s character\n",argv[0]);
	    return 1;
  	}

  	char input=argv[1][0];

	fprintf(stdout,"The complement of '%c' is '%c'.\n",input, complement(input));
}