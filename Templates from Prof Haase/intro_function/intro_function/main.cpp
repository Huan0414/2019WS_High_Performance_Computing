//  Lecture:  8.3.2018
#include "mylib.h"
#include <cassert>                         // assert()
#include <iostream>
using namespace std;




int main()
{
    float varC, varF;               // declaration of variables (not initialize/defined)
    cout << "Hello world!" << endl;
    cout << " Grad Celsius = ";

    cin  >> varC;                    // Input from terminal
    float varK = c2k(varC);          // call function

    cout << " Kelvin : " << varK  << endl;
    cout << "####################\n";

    c2kf(varC, varK, varF);          //  call function

    cout << " Kelvin    : " << varK  << endl;
    cout << " Fahrenheit: " << varF  << endl;

    return 0;
}
