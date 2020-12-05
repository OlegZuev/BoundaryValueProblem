#include <iomanip>
#include <iostream>
#include <string>
#include <fstream>
#include "Functions.h"
using namespace std;

const bool PRINT_TO_FILE = false;

int main()
{
	ofstream fout("../" + to_string(VARIANT) + "_output" + ".txt");
	ostream& out = PRINT_TO_FILE ? fout : cout;

	out << "Variant: " << VARIANT << endl;

	//shooting_method( out);

	finite_difference_scheme_method_obvious(out);
}