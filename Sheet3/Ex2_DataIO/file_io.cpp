#include <iostream>
#include <string>
#include <fstream>
#include <cassert>
#include <vector>
#include <stdexcept>
#include <algorithm>
#include <omp.h>
using namespace std;

#include "file_io.h"

// [Str10, p.364]
void fill_vector(istream& istr, vector<double>& v)
{
    double d=0;
        
    //#pragma omp declare reduction (OurAppend : vector<double> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end())) initializer (omp_priv(omp_orig))
    //#pragma omp parallel shared(istr) reduction(OurAppend:v)
    //{
		while ( istr >> d) v.push_back(d); // Einlesen
	//}
    if (!istr.eof())
    { // Fehlerbehandlung
        cout << " Error handling \n";
        if ( istr.bad() )  throw runtime_error("Schwerer Fehler in istr");
        if ( istr.fail() )   // Versuch des Aufraeumens
        {
            cout << " Failed in reading all data.\n";
            istr.clear();
        }
    }
    v.shrink_to_fit();                 // C++11
        
    return;
}

void read_vector_from_file(const string& file_name, vector<double>& v)
{
    ifstream fin(file_name);           // Oeffne das File im ASCII-Modus
    if( fin.is_open() )                // File gefunden:
    {
       v.clear();                      // Vektor leeren
       fill_vector(fin, v);
    }
    else                               // File nicht gefunden:
    {
        cout << "\nFile " << file_name << " has not been found.\n\n" ;
        assert( fin.is_open() && "File not found." );       // exeption handling for the poor programmer
    }

 return;
}


void write_vector_to_file(const string& file_name, const vector<double>& v)
{
    ofstream fout(file_name);          // Oeffne das File im ASCII-Modus
    if( fout.is_open() )
    {
       for (unsigned int k=0; k<v.size(); ++k)
       {
          fout << v.at(k) << endl;
       }
    }
    else
    {
        cout << "\nFile " << file_name << " has not been opened.\n\n" ;
        assert( fout.is_open() && "File not opened."  );         // exeption handling for the poor programmer
    }

 return;
}

