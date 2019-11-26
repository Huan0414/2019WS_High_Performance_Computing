#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
#include <iterator>
using namespace std;

//int get_primes(unsigned long n);
int single_goldback(long n);
vector<int> count_goldback(long n);

int main(void)
{
    vector <long> primes;
    long n;
    cout << "Please enter a even number larger than 3."<<"\n";
    cin  >> n;
    cout << endl;

    //single number decompositions conditions
    int deNr = single_goldback(n);
    cout << "Overall decompositions number for" <<" "<< n <<" "<<"is: " <<deNr << endl;
    cout << endl;

    //all numbers between 4-n decompositions conditions
    double time1 = 0.0, tstart; //time measurement
    tstart  = clock();
    vector<int> deVector = count_goldback(n);
    time1 = clock() - tstart;
    time1 = time1 /CLOCKS_PER_SEC;

    //Check the maximum
    cout<<"The number has most pairs is" <<" "<< distance(deVector.begin(), max_element(deVector.begin(),deVector.end())+3)<<endl;
    cout<<endl;
    cout<<"Count_goldbach function takes"<< " "<< time1 <<" "<<"Seconds when the range is 4 to" << " "<<n <<endl;

    return 0;
}
