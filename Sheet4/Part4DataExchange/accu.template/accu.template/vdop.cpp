#include "vdop.h"
#include <cassert>               // assert()
#include <cmath>
#include <iostream>
#include <vector>
using namespace std;


void vddiv(vector<double> & x, vector<double> const& y,
                               vector<double> const& z)
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] / z[k];
    }
    return;
}

//******************************************************************************

void vdaxpy(std::vector<double> & x, std::vector<double> const& y,
                       double alpha, std::vector<double> const& z )
{
    assert( x.size()==y.size() && y.size()==z.size() );
    size_t n = x.size();

#pragma omp parallel for
    for (size_t k = 0; k < n; ++k)
    {
        x[k] = y[k] + alpha * z[k];
    }
    return;
}
//******************************************************************************

double dscapr(std::vector<double> const& x, std::vector<double> const& y)
{
    assert( x.size()==y.size());
    size_t n = x.size();

    double    s = 0.0;
//#pragma omp parallel for reduction(+:s)
    for (size_t k = 0; k < n; ++k)
    {
        s += x[k] * y[k];
    }

    return s;
}

//******************************************************************************
//void DebugVector(vector<double> const &v)
//{
    //cout << "\nVector  (nnode = " << v.size() << ")\n";
    //for (size_t j = 0; j < v.size(); ++j)
    //{
        //cout.setf(ios::right, ios::adjustfield);
        //cout << v[j] << "   ";
    //}
    //cout << endl;

    //return;
//}
//******************************************************************************
bool CompareVectors(std::vector<double> const& x, int const n, double const y[], double const eps)
{
    bool bn = (static_cast<int>(x.size())==n);
    if (!bn)
    {
        cout << "#########   Error: " << "number of elements" << endl;
    }
    //bool bv = equal(x.cbegin(),x.cend(),y);
    bool bv = equal(x.cbegin(),x.cend(),y,
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
    );
    if (!bv)
    {
        assert(static_cast<int>(x.size())==n);
        cout << "#########   Error: " << "values" << endl;
    }
    return bn && bv;
}

//******************************************************************************
double par_scalar(vector<double> const &x, vector<double> const &y, MPI_Comm const& icomm)
{
  const double s = dscapr(x,y);
        double sg;
  MPI_Allreduce(&s,&sg,1,MPI_DOUBLE,MPI_SUM,icomm);

  return(sg);
}

//******************************************************************************
void ExchangeAll(vector<double> const &xin, vector<double> &yout, MPI_Comm const &icomm)
{
    int myrank, numprocs,ierr(-1);
    MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
    int const N=xin.size();
    int const sendcount = N/numprocs;                       // equal sized junks
    assert(sendcount*numprocs==N);                          // really all junk sized?
    assert(xin.size()==yout.size());

    auto sendbuf = xin.data();
    auto recvbuf = yout.data();
    ierr = MPI_Alltoall(sendbuf, sendcount, MPI_DOUBLE,
                        recvbuf, sendcount, MPI_DOUBLE, icomm);
    assert(0==ierr);

    return;
}

//******************************************************************************
void ExchangeAllInPlace(vector<double> &xin, MPI_Comm const &icomm)
{
    int myrank, numprocs,ierr(-1);
    MPI_Comm_rank(icomm, &myrank);                          // my MPI-rank
    MPI_Comm_size(icomm, &numprocs);
    int const N=xin.size();
    int const sendcount = N/numprocs;                       // equal sized junks
    assert(sendcount*numprocs==N);                          // really all junk sized?

    auto sendbuf = xin.data();
    ierr = MPI_Alltoall(MPI_IN_PLACE, sendcount, MPI_DOUBLE,
                        sendbuf, sendcount, MPI_DOUBLE, icomm);
    assert(0==ierr);

    return;
}
