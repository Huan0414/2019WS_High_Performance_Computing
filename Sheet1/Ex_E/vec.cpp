#include <iostream>
#include <vector>
#include <algorithm>
#include <ctime>
using namespace std;

double  vecCreate(int n)
 {
    // create a sorted vector and a list container of length n with ascending entries Xk= k+1
     vector<int> contVec;
     contVec.reserve(n);

     for(int i=0;i<n;i++)
     {
         contVec.push_back(i+1);
     }
     //enlarge the container to 2n length and start to measure time
     double t1 = clock();

    //enlarge the container to 2n length
     for(int  i=0;i<n;i++)
     {
     //generate random numbers between 1~n
        int ins;
        ins = rand()%n + 1;
     // insert this into the container and make sure the container is still sorted
        vector<int>::iterator position;
        position = lower_bound(contVec.begin(), contVec.end(), ins); //find the right position
        contVec.insert(position,ins); // insert to the right position
     }

     // what is the overall time does the process take
     double time_span = clock()- t1;
     // check is the new list sorted or not
    if (is_sorted(contVec.cbegin(),contVec.cend()))
    {
        cout<< "The 2n-length vector container is sorted!"<<endl;
    }
    else
    {
        cout<<"Check the vector function!"<<endl;
    }

    cout<<"The time spend to enlarge the vector container to two times length is: \n"
        << time_span/CLOCKS_PER_SEC<<" "<<"seconds"<<endl;;
    cout<<endl;

    return time_span;
 }
