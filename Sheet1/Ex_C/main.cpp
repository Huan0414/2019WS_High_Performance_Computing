#include <iostream>
#include <vector>
#include <string>
#include <chrono>
#include <ctime>
using namespace std;
using namespace std::chrono;

// put this into a head file
long long int sumLoop(int n);
long long int sumFormula(int n);

int main()
{
  cout <<"This exercise is to add all positive integers less or equal n"
       <<"which are multiple of 3 or 5.\n";
  cout <<"Under this purpose, compare run time between for loop and formula calculation."<<endl;
  cout <<endl;

  // enter n
  int n;
  cout <<"please enter a positive integers.\n";
  cin  >> n;
  // check if two functions are correct
  cout <<"The summation value with for-loop is:"<<sumLoop(n)<<endl;
  cout <<"The summation value with a formula is:"<<sumFormula(n)<<endl;
  cout << endl;
  //run each function a thousand time and compare time difference

    // loop function running time
    high_resolution_clock::time_point t1_loop = high_resolution_clock::now();
    for (unsigned i =1;i<=2000;i++)
    {
        sumLoop(n);
    }
    high_resolution_clock::time_point t2_loop = high_resolution_clock::now();
    duration<double> timeSpan_loop = duration_cast<duration<double>>(t2_loop - t1_loop);

    // formula function running time
    high_resolution_clock::time_point t1_formula = high_resolution_clock::now();
    for (unsigned i =1;i<=2000;i++)
    {
        sumFormula(n);
    }
    high_resolution_clock::time_point t2_formula = high_resolution_clock::now();
    duration<double> timeSpan_formula = duration_cast<duration<double>>(t2_formula - t1_formula);

    // check the difference

    if ((timeSpan_loop.count() - timeSpan_formula.count())>0)
    {
        cout<<"Formula function is faster!"<<endl;
        cout<<"The time difference of running 2000 times is: "<< \
        timeSpan_loop.count() - timeSpan_formula.count()<<"seconds."<<endl;
    }
    else
    {
        cout<<"For-loop function is faster!"<<endl;
        cout<<"The time difference of running 2000 times is: "<< \
        timeSpan_formula.count() - timeSpan_loop.count()<<"seconds."<<endl;
    }

    return 0;
}

