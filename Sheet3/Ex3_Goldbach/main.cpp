#include <iostream>
#include <algorithm> 
#include <vector>
#include <ctime>
using namespace std;

#include <mylib.h>

int main()
{
    int n;
    cout << "Please enter a even number larger than 3."<<"\n";
    cin  >> n;
    cout << endl;

    //all numbers between 4-n decompositions conditions
    double t1, tstart; //time measurement
    
	int N = (n-4)/2 + 1;        //Only consider even numbers from 4-n
    vector<int> pair_cont(N,0); // all counts are initialized as zero
    
    tstart  = clock();
    count_goldback(n, pair_cont); 
    
    t1 = clock() - tstart;
    t1 /= CLOCKS_PER_SEC;

    //Check the number with maximum pairs
    cout << "The number has most pairs is" <<" "<< 2*(2 + distance(pair_cont.begin(), max_element(pair_cont.begin(),pair_cont.end())))<< endl;
    cout << endl;
    //Evaluate the performance
  
    

    return 0;
}
