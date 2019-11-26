#include <iostream>
#include <vector>
#include <string>
using namespace std;

long long int sumLoop(int n)

 {
   long long int sumLoop = 0;
    for (int i=1;i<=n;++i)
     {
        if (i%3==0||i%5==0)
        {
        sumLoop += i;
//        cout << sumLoop<<endl;
        }
     }

    return sumLoop;
 }


