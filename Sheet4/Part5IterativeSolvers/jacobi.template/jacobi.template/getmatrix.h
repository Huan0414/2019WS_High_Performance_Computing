#ifndef GETMATRIX_FILE
#define GETMATRIX_FILE

#include "geom.h"
#include <cassert>
#include <vector>
// #####################################################################
/**
 * Abstract matrix class.
 */
class Matrix
{
    public:
       /**
		 * Constructor for abstract matrix class.
         *
         * No memory is allocated.
		 *
		 * @param[in] nrows   number of matrix rows.
         * @param[in] ncols   number of matrix columns.
		*/
       Matrix(int nrows, int ncols);
       //Matrix();

       Matrix(Matrix const &) = default;
       /**
		 * Destructor.
         *
         * No memory is allocated.
		*/
       virtual ~Matrix();

       /**
		 * Checks whether the matrix is a square matrix.
		 *
		 * @return True iff square matrix.
		*/
       bool isSquare() const
       { return _nrows==_ncols;}

       /**
		 * Number of rows in matrix.
		 * @return number of rows.
		 */
       int Nrows() const
          {return _nrows;}

       /**
		 * Number of columns in matrix.
		 * @return number of columns.
		 */
       int Ncols() const
          {return _ncols;}

       /**
		 * Show the matrix entries.
		 */
       virtual void Debug() const = 0;

       /**
        * Extracts the diagonal elements of an inherited matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       virtual void GetDiag(std::vector<double> &d) const = 0;

       /**
        * Extracts the diagonal elements of the matrix.
        *
        * @return  d  vector of diagonal elements
        */
       std::vector<double> const & GetDiag() const
       {
           if ( 0==_dd.size() )          // GH: better?   Nrows()>static_cast<int>(_dd.size())
           {
               _dd.resize(Nrows());
               this->GetDiag(_dd);
            }
            assert( Nrows()==static_cast<int>(_dd.size()) );
            return _dd;
        }

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       virtual void Mult(std::vector<double> &w, std::vector<double> const &u) const = 0;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       virtual void Defect(
                   std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const = 0;

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       virtual int fetch(int row, int col) const =0;

    protected:
       int _nrows;              //!< number of rows in matrix
       int _ncols;              //!< number of columns in matrix
       mutable std::vector<double> _dd; //!< diagonal matrix elements
};

// #####################################################################
class BisectInterpolation;  // class forward declaration
/**
 * Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */
class CRS_Matrix: public Matrix
{
    public:
       /**
        * Constructor
        *
        */
       CRS_Matrix();

       CRS_Matrix(CRS_Matrix const &) = default;


      /**
        * Destructor.
        */
       virtual ~CRS_Matrix() override;
       /**
        * Extracts the diagonal elements of the sparse matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const override;
       ///**
        //* Extracts the diagonal elements of the sparse matrix.
        //*
        //* @return  d  vector of diagonal elements
        //*/
       //std::vector<double> const & GetDiag() const override;

       /**
        * Performs the matrix-vector product  w := K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     u vector
        */
       void Mult(std::vector<double> &w, std::vector<double> const &u) const override;

        /**
        * Calculates the defect/residuum w := f - K*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const override;

       /**
		 * Show the matrix entries.
		 */
       void Debug() const override;

       /**
		 * Finds in a CRS matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const override;

        /**
        * Compare @p this CRS matrix with an external CRS matrix stored in C-Style.
        *
        * The method prints statements on differences found.
        *
        * @param[in]     nnode  row number of external matrix
        * @param[in]     id     start indices of matrix rows of external matrix
        * @param[in]     ik     column indices of external matrix
        * @param[in]     sk     non-zero values of external matrix
        *
        * @return true iff all data are identical.
        */
       bool Compare2Old(int nnode, int const id[], int const ik[], double const sk[]) const;
              
       /**
		 * Calculates the defect and projects it to the next coarser level @f$ f_C := P^T \cdot (f_F - SK\cdot u_F) @f$.  
		 *
		 * @param[in] SK	matrix on fine mesh
		 * @param[in] P	    prolongation operator
		 * @param[in,out] fc  resulting coarse mesh vector (preallocated)
		 * @param[in] ff	r.h.s. on fine mesh
		 * @param[in] uf	status vector on fine mesh 
		 *
		*/
       friend void DefectRestrict(CRS_Matrix const & SK, BisectInterpolation const& P, 
       std::vector<double> &fc, std::vector<double> &ff, std::vector<double> &uf);

    protected:
       //int _nrows;              //!< number of rows in matrix
       //int _ncols;              //!< number of columns in matrix
       int _nnz;                //!< number of non-zero entries
       std::vector<int> _id;    //!< start indices of matrix rows
       std::vector<int> _ik;    //!< column indices
       std::vector<double> _sk; //!< non-zero values
};


/**
 * FEM Matrix in CRS format (compressed row storage; also named CSR),
 * see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 */
class FEM_Matrix: public CRS_Matrix
{
    public:
       /**
        * Initializes the CRS matrix structure from the given discretization in @p mesh.
        *
        * The sparse matrix pattern is generated but the values are 0.
        *
        * @param[in] mesh given discretization
        *
        * @warning A reference to the discretization @p mesh is stored inside this class.
        *          Therefore, changing @p mesh outside requires also
        *          to call method @p Derive_Matrix_Pattern explicitly.
        *
        * @see Derive_Matrix_Pattern
        */
       explicit FEM_Matrix(Mesh const & mesh);

       FEM_Matrix(FEM_Matrix const &) = default;

      /**
        * Destructor.
        */
       ~FEM_Matrix() override;

       /**
        * Generates the sparse matrix pattern and overwrites the existing pattern.
        *
        * The sparse matrix pattern is generated but the values are 0.
       */
       void Derive_Matrix_Pattern()
       {
           //Derive_Matrix_Pattern_slow();
           Derive_Matrix_Pattern_fast();
       }
       void Derive_Matrix_Pattern_fast();
       void Derive_Matrix_Pattern_slow();


        /**
        * Calculates the entries of f.e. stiffness matrix and load/rhs vector @p f for the Laplace operator in 2D.
        * No memory is allocated.
        *
        * @param[in,out] f (preallocated) rhs/load vector
        */
       void CalculateLaplace(std::vector<double> &f);

       /**
        * Applies Dirichlet boundary conditions to stiffness matrix and to load vector @p f.
        * The <a href="https://www.jstor.org/stable/2005611?seq=1#metadata_info_tab_contents">penalty method</a>
        * is used for incorporating the given values @p u.
        *
        * @param[in]     u (global) vector with Dirichlet data
        * @param[in,out] f load vector
        */
       void ApplyDirichletBC(std::vector<double> const &u, std::vector<double> &f);

       ///**
        //* Extracts the diagonal elements of the sparse matrix.
        //*
        //* @param[in,out]  d  (prellocated) vector of diagonal elements
        //*/
       //void GetDiag(std::vector<double> &d) const;   // override in MPI parallel


      /**
        * Adds the element stiffness matrix @p ske and the element load vector @p fe
        * of one triangular element with linear shape functions to the appropriate positions in
        * the stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik).
        *
        * @param[in]     ial   node indices of the three element vertices
        * @param[in]     ske   element stiffness matrix
        * @param[in]     fe    element load vector
        * @param[in,out] f	   distributed local vector storing the right hand side
        *
        * @warning Algorithm assumes  linear triangular elements (ndof_e==3).
       */
       void AddElem_3(int const ial[3], double const ske[3][3], double const fe[3], std::vector<double> &f);


    private:
       Mesh const & _mesh;      //!< reference to discretization

};


///**
 //* Prolongation matrix in CRS format (compressed row storage; also named CSR),
 //* see an <a href="https://en.wikipedia.org/wiki/Sparse_matrix">introduction</a>.
 //*
 //* The prolongation is applied for each node from the coarse mesh to the fine mesh and
 //* is derived only geometrically (no operator weighted prolongation).
 //*/
//class Prolongation: public CRS_Matrix
//{
    //public:
       ///**
        //* Intializes the CRS matrix structure from the given discetization in @p mesh.
        //*
        //* The sparse matrix pattern is generated but the values are 0.
        //*
        //* @param[in] cmesh coarse mesh
        //* @param[in] fmesh fine mesh
        //*
        //* @warning A reference to the discretizations @p fmesh  @p cmesh are stored inside this class.
        //*          Therefore, changing these meshes outside requires also
        //*          to call method @p Derive_Matrix_Pattern explicitely.
        //*
        //* @see Derive_Matrix_Pattern
        //*/
       //Prolongation(Mesh const & cmesh, Mesh const & fmesh);

       ///**
        //* Destructor.
        //*/
       //~Prolongation() override
       //{}

       ///**
        //* Generates the sparse matrix pattern and overwrites the existing pattern.
        //*
        //* The sparse matrix pattern is generated but the values are 0.
       //*/
       //void Derive_Matrix_Pattern() override;

    //private:
       //Mesh const & _cmesh;      //!< reference to coarse discretization
       //Mesh const & _fmesh;      //!< reference to fine discretization
//};

// *********************************************************************


/**
 * Interpolation matrix for prolongation coarse mesh (C) to a fine mesh (F)
 * generated by bisecting edges.
 *
 * All interpolation weights are 0.5 (injection points contribute twice).
 */
class BisectInterpolation: public Matrix
{
    public:
       /**
        * Generates the interpolation matrix for prolongation coarse mesh to a fine mesh
        * generated by bisecting edges.
        * The interpolation weights are all 0.5.
        *
        * @param[in] fathers vector[nnodes][2] containing
        *                    the two coarse grid fathers of a fine grid vertex
        *
        */
       explicit BisectInterpolation(std::vector<int> const & fathers);
       BisectInterpolation();

       BisectInterpolation(BisectInterpolation const &) = default;

       /**
        * Destructor.
        */
       ~BisectInterpolation() override;

       /**
        * Extracts the diagonal elements of the matrix.
        *
        * @param[in,out]  d  (prellocated) vector of diagonal elements
        */
       void GetDiag(std::vector<double> &d) const override;
       ///**
        //* Extracts the diagonal elements of the sparse matrix.
        //*
        //* @return  d  vector of diagonal elements
        //*/
       //std::vector<double> const & GetDiag() const override;

       /**
        * Performs the prolongation  @f$ w_F := P*u_C @f$.
        *
        * @param[in,out] wf resulting fine vector (preallocated)
        * @param[in]     uc coarse vector
        */
       void Mult(std::vector<double> &wf, std::vector<double> const &uc) const override;

       /**
        * Performs the restriction  @f$ u_C := P^T*w_F @f$.
        *
        * @param[in]         wf fine vector
        * @param[in,out]     uc resulting coarse vector (preallocated)
        */
       void MultT(std::vector<double> const &wf, std::vector<double> &uc) const;

        /**
        * Calculates the defect/residuum w := f - P*u.
        *
        * @param[in,out] w resulting vector (preallocated)
        * @param[in]     f load vector
        * @param[in]     u coarse vector
        */
       void Defect(std::vector<double> &w,
                   std::vector<double> const &f, std::vector<double> const &u) const override;

       /**
		 * Show the matrix entries.
		 */
       void Debug() const override;

       /**
		 * Finds in this matrix the access index for an entry at row @p row and column @p col.
		 *
		 * @param[in] row	row index
		 * @param[in] col	column index
		 * @return index for element (@p row, @p col). If no appropriate entry exists then -1 will be returned.
		 *
		 * @warning assert() stops the function in case that matrix element (@p row, @p col) doesn't exist.
		*/
       int fetch(int row, int col) const override;

       /**
		 * Calculates the defect and projects it to the next coarser level @f$ f_C := P^T \cdot (f_F - SK\cdot u_F) @f$.  
		 *
		 * @param[in] SK	matrix on fine mesh
		 * @param[in] P	    prolongation operator
		 * @param[in,out] fc  resulting coarse mesh vector (preallocated)
		 * @param[in] ff	r.h.s. on fine mesh
		 * @param[in] uf	status vector on fine mesh 
		 *
		*/       
       friend void DefectRestrict(CRS_Matrix const & SK, BisectInterpolation const& P, 
       std::vector<double> &fc, std::vector<double> &ff, std::vector<double> &uf);

    protected:
       std::vector<int> _iv;     //!< fathers[nnode][2] of fine grid nodes, double entries denote injection points
       std::vector<double> _vv;  //!< weights[nnode][2] of fathers for grid nodes
};

/**
 * Interpolation matrix for prolongation from coarse mesh (C)) to a fine mesh (F)
 * generated by bisecting edges.
 * 
 * We take into account that values at Dirichlet nodes have to be preserved, i.e.,
 * @f$ w_F = P \cdot I_D \cdot w_C @f$ and @f$ d_C = I_D  \cdot P^T \cdot  d_F@f$
 * with @f$ I_D @f$ as @f$ n_C \times n_C @f$ diagonal matrix and entries
 * @f$ I_{D(j,j)} := \left\{{\begin{array}{l@{\qquad}l} 0 & x_{j}\;\;  \textrm{is Dirichlet node} \\ 1 & \textrm{else} \end{array}}\right. @f$
 *
 * Interpolation weights are eighter 0.5 or 0.0 in case of coarse Dirichlet nodes
 * (injection points contribute twice),
 * Sets weight to zero iff (at least) one father nodes is a Dirichlet node.
 */
class BisectIntDirichlet: public BisectInterpolation
{
    public:
       /**
		 * Default constructor.
		*/
       BisectIntDirichlet()
        : BisectInterpolation()
       {}

       /**
		 * Constructs interpolation from father-@p row and column @p col.
		 *
		 * @param[in] fathers	two father nodes from each fine node [nnode_f*2].
		 * @param[in] idxc_dir	vector containing the indices of coarse mesh Dirichlet nodes.
		 *
		*/
       BisectIntDirichlet(std::vector<int> const & fathers, std::vector<int> const & idxc_dir);

       BisectIntDirichlet(BisectIntDirichlet const &) = default;


       /**
        * Destructor.
        */
       ~BisectIntDirichlet() override;
};



// *********************************************************************

/**
 * Calculates the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions.
 * @param[in]	ial	node indices of the three element vertices
 * @param[in]	xc	vector of node coordinates with x(2*k,2*k+1) as coordinates of node k
 * @param[out] ske	element stiffness matrix
 * @param[out] fe	element load vector
 */
void CalcElem(int const ial[3], double const xc[], double ske[3][3], double fe[3]);

/**
 * Adds the element stiffness matrix @p ske and the element load vector @p fe
 * of one triangular element with linear shape functions to the appropriate positions in
 * the symmetric stiffness matrix, stored as CSR matrix K(@p sk,@p id, @p ik)
 *
 * @param[in] ial   node indices of the three element vertices
 * @param[in] ske	element stiffness matrix
 * @param[in] fe	element load vector
 * @param[out] sk	vector non-zero entries of CSR matrix
 * @param[in] id	index vector containing the first entry in a CSR row
 * @param[in] ik	column index vector of CSR matrix
 * @param[out] f	distributed local vector storing the right hand side
 *
 * @warning Algorithm requires indices in connectivity @p ial in ascending order.
 *          Currently deprecated.
*/
void AddElem(int const ial[3], double const ske[3][3], double const fe[3],
             int const id[], int const ik[], double sk[], double f[]);



#endif
