// see:   http://llvm.org/docs/CodingStandards.html#include-style
#include "vdop.h"
//#include "geom.h"
#include "par_geom.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <ctime>                  // contains clock()
#include <fstream>
#include <iostream>
#include <list>
#include <numeric>                // accumulate()
#include <string>
#include <vector>

using namespace std;


ParMesh::ParMesh(int ndim, int nvert_e, int ndof_e, int nedge_e, MPI_Comm const &icomm)
    : Mesh(ndim, nvert_e, ndof_e, nedge_e),
      _icomm(icomm), _numprocs(-1), _myrank(-1),
      _v_l2g(0), _t_l2g(0), _v_g2l{{}}, _t_g2l{{}}, _valence(0),
      _sendbuf(0), _sendcounts(0), _sdispls(0),
      _loc_itf(0), _gloc_itf(0), _buf2loc(0)
{
    MPI_Comm_size(icomm, &_numprocs);
    MPI_Comm_rank(icomm, &_myrank);
}

ParMesh::~ParMesh()
{}



ParMesh::ParMesh(string const &sname, MPI_Comm const &icomm)
    : ParMesh(2, 3, 3, 3, icomm) // two dimensions, 3 vertices, 3 dofs, 3 edges per element
{
    //const int numprocs = _icomm.Get_size();
    const string NS    = "_" + to_string(_numprocs);
    const string fname = sname + NS + ".txt";
    //cout << "############  " << fname << endl;
    ReadVertexBasedMesh(fname);
    // ------------------------------------------------------------------------------
    // Until this point  a l l  processes possess  a l l  mesh info in  g l o b a l  numbering
    //
    // Now, we have to select the data belonging to my_rank
    // and we have to create the mapping local to global (l2g) and vice versa (g2l)
    // ------------------------------------------------------------------------------

    // save the global node mesh (maybe we need it later)
    DeriveEdgeFromVertexBased();                     //                    and even more
    Mesh global_mesh(*this);                         // requires a  l o t  of memory
    Del_EdgeConnectivity();

    // read the subdomain info
    const string dname = sname + NS + "_sd" + ".txt";
    vector<int> t2d = ReadElementSubdomains(dname);  // global mapping triangle to subdomain for all elements

    //const int myrank   = _icomm.Get_rank();
    Transform_Local2Global_Vertex(_myrank, t2d);      // Vertex based mesh: now in  l o c a l  indexing

    DeriveEdgeFromVertexBased();                     // Generate also the  l o c a l  edge based information

    Generate_VectorAdd();


    // Now we have to organize the MPI communication of vertices on the subdomain interfaces

    return;
}

vector<int> ParMesh::ReadElementSubdomains(string const &dname)
{
    ifstream ifs(dname);
    if (!(ifs.is_open() && ifs.good())) {
        cerr << "ParMesh::ReadElementSubdomain: Error cannot open file " << dname << endl;
        assert(ifs.is_open());
    }

    int const OFFSET(1);             // Matlab to C indexing
    cout << "ASCI file  " << dname << "  opened" << endl;

    // Read some mesh constants
    int nelem;
    ifs >> nelem;
    cout << nelem << endl;
    assert( Nelems() == nelem);

    // Allocate memory
    vector<int> t2d(nelem, -1);
    // Read element mapping
    for (int k = 0; k < nelem; ++k) {
        int tmp;
        ifs >> tmp;
        t2d[k] = tmp - OFFSET;
    }

    return t2d;
}

void ParMesh::Transform_Local2Global_Vertex(int const myrank, vector<int> const &t2d)
{
    // number of local elements
    const int l_ne = count(t2d.cbegin(), t2d.cend(), myrank);
    //cout << myrank <<  "::  " << lne << endl;
    vector<int> l_ia(l_ne * NverticesElements(), -1); // local elements still with global vertex numbers
    _t_l2g.resize(l_ne, -1);

    int lk = 0;
    for (size_t k = 0; k < t2d.size(); ++k) {
        if (myrank == t2d[k]) {
            //if (0==myrank)
            //{
            //cout << lk << "   k  " <<  t2d[k] << endl;
            //}
            l_ia[3 * lk  ] = _ia[3 * k  ];
            l_ia[3 * lk + 1] = _ia[3 * k + 1];
            l_ia[3 * lk + 2] = _ia[3 * k + 2]; // local elements still with global vertex numbers
            //l_ia.at(3*lk+2) = _ia.at(3*k+2);
            _t_l2g[lk] = k;                       // elements: local to global mapping
            _t_g2l[k] = lk;                       //           global to local
            ++lk;
        }
    }
    // Checks:
    assert( count(l_ia.cbegin(), l_ia.cend(), -1)   == 0 );
    assert( count(_t_l2g.cbegin(), _t_l2g.cend(), -1) == 0 );

    // Vertices: local to global mapping
    auto tmp = l_ia;
    sort(tmp.begin(), tmp.end());
    auto ip = unique(tmp.begin(), tmp.end());
    tmp.erase(ip, tmp.end());
    _v_l2g = tmp;                                 // Vertices: local to global mapping
    for (size_t lkv = 0; lkv < _v_l2g.size(); ++lkv) {
        _v_g2l[_v_l2g[lkv]] = lkv;                //           global to local
    }

    // Boundary edges
    vector<int> l_bedges;
    for (size_t b = 0; b < _bedges.size(); b += 2) {
        int const v1 = _bedges[b];                // global vertex numbers
        int const v2 = _bedges[b + 1];
        try {
            int const lv1 = _v_g2l.at(v1);         // map[] would add that element
            int const lv2 = _v_g2l.at(v2);         //       but at() throws an exeption
            l_bedges.push_back(lv1);
            l_bedges.push_back(lv2);              // Boundaries: already in local indexing
        }
        catch (std::out_of_range & err) {
            //cout << ".";
        }
    }

    // number of local vertices
    const int l_nn = _v_l2g.size();
    vector<double> l_xc(Ndims()*l_nn);
    for (int lkk = 0; lkk < l_nn; ++lkk) {
        int k = _v_l2g.at(lkk);
        l_xc[2 * lkk  ]   = _xc[2 * k  ];
        l_xc[2 * lkk + 1] = _xc[2 * k + 1];
    }


    // Now, we represent the vertex mesh in  l o c a l  numbering
    // elements

    for (size_t i = 0; i < l_ia.size(); ++i) {
        l_ia[i] = _v_g2l.at(l_ia[i]);              // element vertices: global to local
    }
    SetNelem(l_ne);
    _ia = l_ia;
    // boundary
    _bedges = l_bedges;
    // coordinates
    SetNnode(l_nn);
    _xc = l_xc;

    return;
}


void ParMesh::Generate_VectorAdd()
{
    // Some checks
    int lnn = Nnodes();                           // local number of vertices
    assert(static_cast<int>(_v_l2g.size()) == lnn);
    int ierr;

    // ---- Determine global largest vertex index
    int gidx_max{-1};                             // global largest vertex index
    int lmax = *max_element(_v_l2g.cbegin(), _v_l2g.cend());
    MPI_Allreduce(&lmax, &gidx_max, 1, MPI_INT, MPI_MAX, _icomm);
    int gidx_min{-1};                             // global smallest vertex index
    int lmin = *min_element(_v_l2g.cbegin(), _v_l2g.cend());
    MPI_Allreduce(&lmin, &gidx_min, 1, MPI_INT, MPI_MIN, _icomm);
    //cout << gidx_min << "  " << gidx_max << endl;
    assert(0 == gidx_min);                        // global indices have to start with 0


    // ---- Determine for all global vertices the number of subdomains it belongs to
    vector<int> global(gidx_max, 0);              // global scalar array for vertices
    for (auto const gidx : _v_l2g)  global[gidx] = 1;
    // https://www.mpi-forum.org/docs/mpi-2.2/mpi22-report/node109.htm
    ierr = MPI_Allreduce(MPI_IN_PLACE, global.data(), global.size(), MPI_INT, MPI_SUM, _icomm);
    if (0 == MyRank())  cout << global << endl;
    //MPI_Barrier(_icomm);
    //cout << _xc[2*_v_g2l.at(2)] << " , " << _xc[2*_v_g2l.at(2)+1] << endl;
    //MPI_Barrier(_icomm);

    //  now, global[] contains the number of subdomains a global vertex belongs to
    if ( count(global.cbegin(), global.cend(), 0) > 0 )
        cerr << "\n !!!   Non-continuous global vertex indexing  !!!\n";

    // ---- Determine local interface vertices ( <==> global[] > 1 )
    //      _loc_itf,  neigh_itf
    //vector<int> loc_itf;                          // local indices of interface vertices on this MPI process
    for (size_t lk = 0; lk < _v_l2g.size(); ++lk) {
        int const gk = _v_l2g[lk];                // global index of local vertex lk
        if ( global[gk] > 1 ) {
            _loc_itf.push_back(lk);               // local indices of interface vertices on this MPI process
        }
    }

    MPI_Barrier(_icomm);
    if (0 == MyRank())  cout << "\n..._loc_itf...\n" << _loc_itf << "\n......\n";
    MPI_Barrier(_icomm);
    // ---- global indices of local interface vertices
    //auto gloc_itf(_loc_itf);
    _gloc_itf=_loc_itf;
    for_each(_gloc_itf.begin(), _gloc_itf.end(), [this] (auto & v) -> void { v = _v_l2g[v];} );
    MPI_Barrier(_icomm);
    if (0 == MyRank())  cout << "\n..._gloc_itf...\n" << _gloc_itf << "\n......\n";
    //DebugVector(_gloc_itf,"_gloc_itf");

    // ---- Determine the global length of interfaces
    vector<int> vnn(NumProcs(), -1);              // number of interface vertices per MPI rank
    int l_itf(_loc_itf.size());                    // # local interface vertices
    ierr = MPI_Allgather(&l_itf, 1, MPI_INT, vnn.data(), 1, MPI_INT, _icomm);
    assert(0 == ierr);
    //cout << vnn << endl;

    // ---- Now we consider only the inferface vertices
    int snn = accumulate(vnn.cbegin(), vnn.cend(), 0); // required length of array for global interface indices
    //cout << snn << "   " << gnn << endl;
    vector<int> dispnn(NumProcs(), 0) ;           // displacement of interface vertices per MPI rank
    partial_sum(vnn.cbegin(), vnn.cend() - 1, dispnn.begin() + 1);
    //cout << dispnn << endl;

    // ---- Get the global indices for all global interfaces
    vector<int> g_itf(snn, -1);                    // collects all global indices of the global interfaces
    // https://www.mpich.org/static//docs/v3.0.x/www3/MPI_Gatherv.html
    ierr = MPI_Gatherv( _gloc_itf.data(), _gloc_itf.size(), MPI_INT,
                        g_itf.data(), vnn.data(), dispnn.data(), MPI_INT, 0, _icomm);
    assert(0 == ierr);
    // https://www.mpich.org/static/docs/v3.1/www3/MPI_Bcast.html
    ierr = MPI_Bcast(g_itf.data(), g_itf.size(), MPI_INT, 0, _icomm);
    assert(0 == ierr);                            // Now, each MPI rank has the all global indices of the global interfaces
    MPI_Barrier(_icomm);
    if (MyRank() == 0) cout << "\n...g_itf...\n" << g_itf << "\n......\n";
    MPI_Barrier(_icomm);

    // ----- Determine all MPI ranks a local interface vertex belongs to
    vector<vector<int>> neigh_itf(_loc_itf.size());// subdomains a local interface vertex belongs to
    for (size_t lk = 0; lk < _loc_itf.size(); ++lk) {
        const int gvert = _gloc_itf[lk];           // global index of local interface node lk
        for (int rank = 0; rank < NumProcs(); ++rank) {
            auto const startl = g_itf.cbegin() + dispnn[rank];
            auto const endl   = startl + vnn[rank];
            if ( find( startl, endl, gvert) != endl) {
                neigh_itf[lk].push_back(rank);
            }
        }
    }

    // ---- check the available info in  _loc_itf[lk], _gloc_itf[lk], neigh_itf[lk]
    MPI_Barrier(_icomm);
    //if (MyRank()==0)  cout << "\n...neigh_itf ...\n"  << neigh_itf << endl;
    if (MyRank() == 0) {
        for (size_t lk = 0; lk < _loc_itf.size(); ++lk ) {
            cout << lk << " : local idx " << _loc_itf[lk] << " , global idx " <<  _gloc_itf[lk];
            cout << "  with MPI ranks " << neigh_itf[lk] << endl;
        }
    }
    MPI_Barrier(_icomm);
    
    // ---- store the valence (e.g., the number of subdomains it belongs to) of all local vertices
    _valence.resize(Nnodes(),1);
    for (size_t lk = 0; lk < _loc_itf.size(); ++lk)
    {
        _valence[_loc_itf[lk]] = neigh_itf[lk].size();
    }
    //DebugVector(_valence,"_valence",_icomm);

    // ---- We ware going to use MPI_Alltoallv for data exchange on interfaces
    //      https://www.mpi-forum.org/docs/mpi-3.1/mpi31-report/node109.htm#Node109
    //      https://www.open-mpi.org/doc/v4.0/man3/MPI_Alltoallv.3.php
    //int MPI_Alltoallv(const void* sendbuf, const int sendcounts[], const int sdispls[], MPI_Datatype sendtype, void* recvbuf, const int recvcounts[], const int rdispls[], MPI_Datatype recvtype, MPI_Comm comm)
    //
    //  MPI_Alltoallv needs:
    //      vector<double> sendbuf             (MPI_IN_PLACE: used also as recvbuf)
    //      vector<int> sendcounts             (the same as for recv)
    //      vector<int> sdispls                (the same as for recv)
    //
    //  We need to map the interface vertices onto the sendbuffer:
    //      vector<int> loc_itf                 local  index of interface vertex lk
    //      vector<int> gloc_itf                global index of interface vertex lk
    //      vector<int> buf2loc                 local indices of sendbuffer positions (the same as for recv)

    // ---- Determine sendcounts[] and sdipls[] from neigh_itf[]
    //vector<int> _sendcounts(NumProcs(), 0);
    _sendcounts.resize(NumProcs(), 0);
    for (size_t lk = 0; lk < _loc_itf.size(); ++lk ) {
        auto const &kneigh = neigh_itf[lk];
        for (size_t ns = 0; ns < kneigh.size(); ++ns) {
            ++_sendcounts[kneigh[ns]];
        }
    }
    if (MyRank() == 0)  cout << "\n..._sendcounts ...\n"  << _sendcounts << endl;

    //vector<int> _sdispls(NumProcs(), 0);
    _sdispls.resize(NumProcs(), 0);
    partial_sum(_sendcounts.cbegin(), _sendcounts.cend() - 1, _sdispls.begin() + 1);
    //vector<int> _sdispls(NumProcs()+1, 0);
    //partial_sum(_sendcounts.cbegin(), _sendcounts.cend(), _sdispls.begin() + 1);
    if (MyRank() == 0)  cout << "\n..._sdispls ...\n"  << _sdispls << endl;

    // ---- Determine size of buffer 'nbuffer' and mapping 'buf2loc'
    int const nbuffer = accumulate(_sendcounts.cbegin(), _sendcounts.cend(), 0);
    //vector<int> _buf2loc(nbuffer, -1);
    _buf2loc.resize(nbuffer, -1);
    int buf_idx = 0;                              // position in buffer
    for (int rank = 0; rank < NumProcs(); ++rank) {
        assert( buf_idx ==  _sdispls[rank]);
        for (size_t lk = 0; lk < _loc_itf.size(); ++lk ) {
            auto const &kneigh = neigh_itf[lk];
            if (find(kneigh.cbegin(),kneigh.cend(),rank)!=kneigh.cend())
            {
               _buf2loc[buf_idx] =  _loc_itf[lk];
               ++buf_idx;
            }
        }
    }
    //if (MyRank() == 0)  cout << "\n...buf2loc ...\n"  << buf2loc << endl;
    //DebugVector(buf2loc,"buf2loc",_icomm);

    // ---- Allocate send/recv buffer
    //vector<double> _sendbuf(nbuffer,-1.0);
    _sendbuf.resize(nbuffer,-1.0);

    assert(CheckInterfaceExchange_InPlace());
    cout << " Check of data exchange (InPlace) successful!\n";
    assert(CheckInterfaceExchange());
    cout << " Check of data exchange successful!\n";
    assert(CheckInterfaceAdd_InPlace());
    cout << " Check of data add successful!\n";
    assert(CheckInterfaceAdd());
    cout << " Check of data add (InPlace) successful!\n";
    
    vector<double> x(Nnodes(),-1.0);
    VecAccu(x);
    cout << " VecAccu (InPlace) successful!\n";

    
    return;
}

bool ParMesh::CheckInterfaceExchange_InPlace() const
{
    vector<double> x(Nnodes(),-1.0);
    copy(_v_l2g.cbegin(),_v_l2g.cend(),x.begin());       // init x with global vertex indices

    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        _sendbuf[ls] = x[_buf2loc.at(ls)];
    }
    int ierr = MPI_Alltoallv(MPI_IN_PLACE, _sendcounts.data(), _sdispls.data(), MPI_DOUBLE,
                          _sendbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE, _icomm);
    assert(ierr==0);
    //DebugVector(_sendbuf,"_sendbuf",_icomm);

    vector<double> y(x);
    for(size_t lk = 0; lk<_loc_itf.size(); ++lk) y[_loc_itf.at(lk)] = -1.0;   // only for interface nodes
    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        y[_buf2loc.at(ls)] = _sendbuf[ls];
    }

    double const  eps=1e-10;
    bool bv = equal(x.cbegin(),x.cend(),y.cbegin(),
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
                   );
    return bv;
}

bool ParMesh::CheckInterfaceExchange() const
{
    vector<double> x(Nnodes(),-1.0);
    copy(_v_l2g.cbegin(),_v_l2g.cend(),x.begin());       // init x with global vertex indices

    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        _sendbuf[ls] = x[_buf2loc.at(ls)];
    }
    vector<double> recvbuf(_sendbuf.size());
    int ierr = MPI_Alltoallv(_sendbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE,
                              recvbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE, _icomm);
    //DebugVector(_sendbuf,"_sendbuf",_icomm);
    //DebugVector(recvbuf,"recvbuf",_icomm);
    assert(ierr==0);

    vector<double> y(x);
    for(size_t lk = 0; lk<_loc_itf.size(); ++lk) y[_loc_itf.at(lk)] = -1.0;   // only for interface nodes
    for(size_t ls = 0; ls<recvbuf.size(); ++ls)
    {
        y[_buf2loc.at(ls)] = recvbuf[ls];
    }
    //cout << "WRONG  : " << count(y.cbegin(),y.cend(), -1.0) << endl;

    double const  eps=1e-10;
    bool bv = equal(x.cbegin(),x.cend(),y.cbegin(),
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
                   );
    return bv;
}

bool ParMesh::CheckInterfaceAdd_InPlace() const
{
    vector<double> x(Nnodes(),-1.0);
    for (size_t i=0; i<x.size(); ++i)
    {
        x[i] = _xc[2*i]+_xc[2*i+1];                      // init x with coordinate values
    }

    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        _sendbuf[ls] = x[_buf2loc.at(ls)];
    }
    int ierr = MPI_Alltoallv(MPI_IN_PLACE, _sendcounts.data(), _sdispls.data(), MPI_DOUBLE,
                          _sendbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE, _icomm);
    assert(ierr==0);
    //DebugVector(_sendbuf,"_sendbuf",_icomm);

    vector<double> y(x);
    for(size_t lk = 0; lk<_loc_itf.size(); ++lk) y[_loc_itf.at(lk)] = 0.0;   // only for interface nodes
    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        y[_buf2loc.at(ls)] += _sendbuf[ls];
    }
    MPI_Barrier(_icomm);
    //DebugVector(x,"x",_icomm);
    //DebugVector(y,"y",_icomm);
    for (size_t i= 0; i<y.size(); ++i) y[i]/=_valence[i];   // divide by valence

    double const  eps=1e-10;
    bool bv = equal(x.cbegin(),x.cend(),y.cbegin(),
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
                   );
    return bv;
}

bool ParMesh::CheckInterfaceAdd() const
{
    vector<double> x(Nnodes(),-1.0);
    for (size_t i=0; i<x.size(); ++i)
    {
        //x[i] = _xc[2*i]+_xc[2*i+1];                      // init x with coordinate values
        x[i] = _v_l2g[i];
    }
    
    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        _sendbuf[ls] = x[_buf2loc.at(ls)];
    }
    vector<double> recvbuf(_sendbuf.size());
    int ierr = MPI_Alltoallv(_sendbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE,
                              recvbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE, _icomm);
    //DebugVector(_sendbuf,"_sendbuf",_icomm);
    //DebugVector(recvbuf,"recvbuf",_icomm);
    assert(ierr==0);
 
    vector<double> y(x);
    for(size_t lk = 0; lk<_loc_itf.size(); ++lk) y[_loc_itf.at(lk)] = 0.0;   // only for interface nodes
    for(size_t ls = 0; ls<recvbuf.size(); ++ls)
    {
        //if (0==MyRank()) cout << ls << ": " << _buf2loc.at(ls) << "  " << y[_buf2loc.at(ls)] << "("<< x[_buf2loc.at(ls)] << ")" << "  " << recvbuf[ls] << "  (" <<  _sendbuf[ls] << ")" << endl;
        y[_buf2loc.at(ls)] += recvbuf[ls];
    }
    MPI_Barrier(_icomm);
    //DebugVector(x,"x",_icomm);
    //DebugVector(y,"y",_icomm);
    for (size_t i= 0; i<y.size(); ++i) y[i]/=_valence[i];   // divide by valence

    double const  eps=1e-10;
    bool bv = equal(x.cbegin(),x.cend(),y.cbegin(),
                          [eps](double a, double b) -> bool
                          { return std::abs(a-b)<eps*(1.0+0.5*(std::abs(a)+ std::abs(a))); }
                   );
    return bv;
}


// ----------

void ParMesh::VecAccu(std::vector<double> &w) const
{
    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        _sendbuf[ls] = w[_buf2loc.at(ls)];
    }
    int ierr = MPI_Alltoallv(MPI_IN_PLACE, _sendcounts.data(), _sdispls.data(), MPI_DOUBLE,
                          _sendbuf.data(), _sendcounts.data(), _sdispls.data(), MPI_DOUBLE, _icomm);
    assert(ierr==0);
    //DebugVector(_sendbuf,"_sendbuf",_icomm);

    for(size_t lk = 0; lk<_loc_itf.size(); ++lk) w[_loc_itf.at(lk)] = 0.0;   // only for interface nodes
    for(size_t ls = 0; ls<_sendbuf.size(); ++ls)
    {
        w[_buf2loc.at(ls)] += _sendbuf[ls];
    }

    return;
}







































