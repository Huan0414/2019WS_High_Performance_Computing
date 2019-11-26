#include <iostream>
#include <list>
#include <algorithm>
#include <ctime>

using namespace std;

double  listCreate(int n)
{
    // create a sorted vector and a list container of length n with ascending entries Xk= k+1
     list<int> contList(n);
     for(int i=0;i<n;i++)
     {
         contList.push_back(i+1); // or use iota, should use iterator here
     }

     //enlarge the container to 2n length and start to measure time
     double t1 = clock(); //starter time

     for(int i=0;i<n;i++)
     {
     //generate random numbers between 1~n
        int ins= rand()%n + 1;

     // insert this into the container and make sure the container is still sorted
        list<int>::iterator position = lower_bound(contList.begin(), contList.end(), ins); //find the right position
        contList.insert(position,ins); // insert to the right position
     }

     // what is the overall time does the process take
     double time_span = clock()-t1;

     // check is the new list sorted or not
     if (is_sorted(contList.begin(),contList.end()))
    {
        cout<< "The 2n-length list container is sorted!"<<endl;
    }
    else
    {
        cout<<"Check the list function!"<<endl;
    }

    cout<<"The time spend to enlarge the list container to two times length is: \n"
        << time_span/CLOCKS_PER_SEC<<" "<<"seconds"<<endl;
    cout<<endl;

    return time_span;

}

