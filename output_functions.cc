#include "shell.h"
using namespace blitz;


void writeRealData (Array<double,1> *A, ofstream& file)
{
	for (int i = 2; i < N+2; ++i)
	{
		file<<setw(10)<<(*A)(i)<<"\t\t";
	}
	file<<"\n";
}
