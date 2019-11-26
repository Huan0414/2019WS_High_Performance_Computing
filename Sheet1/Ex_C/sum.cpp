#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
using namespace std;

vector<double> calMean(vector<double> const &x); //put this into a head file

int main()
{
    cout << "File einlesen." << endl;

    const string name("data_1.txt");              // name of input file
    const string name2("out_1.txt");               // name of output file
    vector<double> a;

    read_vector_from_file(name, a);

    vector<double> b(6);
    int N = a.size();
    // determine minimum and maximum
    b.at(0) = *min_element(a.begin(),a.end()); //minimum
    b.at(1) = *max_element(a.begin(),a.end()); //maximum
    vector<double> dMeanArray = calMean(a);
    b.at(2) = dMeanArray[0]; //arithmetic mean
    b.at(3) = dMeanArray[1]; //geometric mean
    b.at(4) = dMeanArray[2]; //harmonic mean
    b.at(5) = dMeanArray[3]; //standard deviation
    // print out
    for (unsigned int k=0; k<b.size(); ++k)
    {
        cout << "  " << b.at(k);
    }
    cout << endl;

    write_vector_to_file(name2, b);

    return 0;
}
