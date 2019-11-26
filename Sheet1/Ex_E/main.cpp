#include <iostream>
#include <list>
#include <vector>
#include <random>
#include <algorithm>
#include <ctime>
#include <chrono>
using namespace std;

double  listCreate(int n);
double  vecCreate(int n);

int main()
{
    cout <<"This exercise is to compare vectors and list"<<endl;
    cout << endl;

    //get n
    int n;
    cout <<"Please enter a number n larger than 1000"<<endl;
    cin  >> n;
    cout << endl;


    if (listCreate(n)>vecCreate(n))
    {
        cout<<"Vector container is faster!\n"
            <<"Because list lacks direct access to the elements by iterate!"<<endl;
    }
    else
    {
        cout <<"Something is wrong!";
    }
    return 0;

}
