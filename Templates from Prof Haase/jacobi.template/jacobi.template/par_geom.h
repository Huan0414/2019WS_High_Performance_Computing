#ifndef PAR_GEOM_FILE
#define PAR_GEOM_FILE
#include "geom.h"
#include "vdop.h"
#include <array>
#include <functional>             // function; C++11
#include <iostream>
#include <map>
#include <memory>                 // shared_ptr
#include <mpi.h>                  // MPI
#include <string>
#include <vector>

class ParMesh: public Mesh
{
        public:
       /**
		 * Constructor initializing the members with default values.
         *
		 * @param[in] ndim  space dimensions (dimension for coordinates)
         * @param[in] nvert_e  number of vertices per element (dimension for connectivity)
         * @param[in] ndof_e   degrees of freedom per element (= @p nvert_e for linear elements)
         * @param[in] nedge_e  number of edges per element (= @p nvert_e for linear elements in 2D)
		 */
       explicit ParMesh(int ndim, int nvert_e=0, int ndof_e=0, int nedge_e=0, MPI_Comm const& icomm=MPI_COMM_WORLD);

       ParMesh(ParMesh const&) = default;

       ParMesh& operator=(ParMesh const&) = delete;


		/**
		 * Destructor.
         *
         * See clang warning on
         * <a href="https://stackoverflow.com/questions/28786473/clang-no-out-of-line-virtual-method-definitions-pure-abstract-c-class/40550578">weak-vtables</a>.
		 */
       virtual ~ParMesh();

		/**
		 * Reads mesh data from a binary file.
         *
		 * @param[in] sname	suffix of file name
         * @param[in] icomm  MPI communicator
         * @see  ascii_write_mesh.m for the file format.
		*/
	   explicit ParMesh(std::string const &sname, MPI_Comm const& icomm=MPI_COMM_WORLD);

       void VecAccu(std::vector<double> &w) const;

       /** 	Inner product
	    * @param[in] x	vector
	    * @param[in] y	vector
	    * @return 	    resulting Euclidian inner product <x,y>
       */
       double dscapr(std::vector<double> const &x, std::vector<double> const &y) const
       {return par_scalar(x, y, _icomm);}

    private:
        /**
		 * Reads the global triangle to subdomain mapping.
         *
		 * @param[in] dname	file name
         *
         * @see ascii_write_subdomains.m for the file format
		*/
       std::vector<int> ReadElementSubdomains(std::string const &dname);


       /**
		 * Transform
         *
		 * @param[in] myrank	MPI rank of this process
         * @param[in] t2d       global mapping triangle to subdomain for all elements (vertex based)
		*/
	   void Transform_Local2Global_Vertex(int myrank, std::vector<int> const &t2d);


       /**
		 * Transform
		*/
       void Generate_VectorAdd();

       bool CheckInterfaceExchange_InPlace() const;
       bool CheckInterfaceExchange() const;
       bool CheckInterfaceAdd_InPlace() const;
       bool CheckInterfaceAdd() const;


    public:
        /** 	MPI rank of the calling process in communication group.
        *
	    * @return 	    MPI rank of the calling process
       */
       int MyRank() const {return _myrank;}

       /** 	Number of MPI processes in communication group.
        *
	    * @return 	    Number of MPI processes
       */
       int NumProcs() const {return _numprocs;}


    private:
       // Don't use  &_icomm  ==> Error
       MPI_Comm const _icomm;                     //!< MPI communicator for the group of processes
       int _numprocs;                             //!< number of MPI processes
       int _myrank;                               //!< my MPI rank
       std::vector<int> _v_l2g;                   //!< vertices:  local to global mapping
       std::vector<int> _t_l2g;                   //!< triangles: local to global mapping
       std::map<int,int> _v_g2l;                  //!< vertices:  global to local mapping
       std::map<int,int> _t_g2l;                  //!< triangles: global to local mapping

       //std::vector<int> e_l2g;                    //!< edges:     local to global mapping

       std::vector<int> _valence;                 //!< valence of local vertices, i.e. number of subdomains they belong to
       //  MPI_Alltoallv needs:
       mutable std::vector<double> _sendbuf;      //!< send buffer  a n d  receiving buffer (MPI_IN_PLACE)
       std::vector<int>    _sendcounts;           //!< number of data to send to each MPI rank (the same as for recv)
       std::vector<int>    _sdispls;              //!< offset of data to send to each MPI rank wrt. _senbuffer (the same as for recv)
    //
    //  We need to map the interface vertices onto the sendbuffer:
       std::vector<int> _loc_itf;                 //!< local  index of interface vertex lk
       std::vector<int> _gloc_itf;                //!< global index of interface vertex lk
       std::vector<int> _buf2loc;                 //!< local indices of sendbuffer positions (the same as for recv)


};


#endif
