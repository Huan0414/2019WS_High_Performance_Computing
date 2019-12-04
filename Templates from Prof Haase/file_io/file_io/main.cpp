#include <iostream>
#include <vector>
#include <string>
using namespace std;

#include "file_io.h"

int main()
{
    cout << "File einlesen." << endl;

    const string name("f_1_d.txt");              // name of input file
    const string name2("out.txt");               // name of output file
    vector<double> a;

    read_vector_from_file(name, a);

    for (unsigned int k=0; k<a.size(); ++k)
    {
        cout << "  " << a.at(k);
    }
    cout << endl;

    write_vector_to_file(name2, a);


    return 0;
}
