#ifndef JACSOLVE_FILE
#define JACSOLVE_FILE
#include "geom.h"
#include "getmatrix.h"
#include <vector>
#include <mpi.h>            // MPI
#include "par_geom.h"


/**
 * Solves linear system of equations  K @p u = @p f  via the Jacobi iteration.
 * We use a distributed symmetric  CSR matrix @p SK and initial guess of the
 * solution is set to 0.
 * @param[in] SK	CSR matrix
 * @param[in] f		distributed local vector storing the right hand side
 * @param[out] u	accumulated local vector storing the solution.
*/
void JacobiSolve(CRS_Matrix const &SK, std::vector<double> const &f, std::vector<double> &u);


void JacobiSolveParallel(ParMesh const &cmesh, CRS_Matrix const &SK, std::vector<double> const &fv, std::vector<double> &uv, MPI_Comm const &icomm);

/**
 * Solves linear system of equations  K @p u = @p f  via the Jacobi iteration.
 * We use a distributed symmetric  CSR matrix @p SK and initial guess of the
 * solution is set to 0.
 *
 * In each smoothing step: @f$ \widehat{u} := u + \omega D^{-1}\left({f-K*u}\right) @f$
 *
 * @param[in]  SK	CSR matrix
 * @param[in]  f	distributed local vector storing the right hand side
 * @param[out] u	accumulated local vector storing the solution.
 * @param[in,out] r	auxiliary local vector.
 * @param[in] nsmooth	number of smoothing steps.
 * @param[in] omega	relaxation parameter.
 * @param[in] zero	initial solution @p u is assumed to be zero.
*/
void JacobiSmoother(Matrix const &SK, std::vector<double> const &f, std::vector<double> &u,
                    std::vector<double> & r, int nsmooth=1, double const omega=1.0, bool zero=false);

/**
 * @brief Simple diagonale preconditioning.
 *
 * The residuum @p r scaled by the inverse diagonal of matríx @p SK results in the correction @p w.
 *
 * @f$ w :=  \omega D^{-1}*r @f$
 *
 * @param[in]  SK	matrix
 * @param[in]  r	distributed local vector storing the residuum
 * @param[out] w	accumulated local vector storing the correction.
 * @param[in] omega	relaxation parameter.
*/
void DiagPrecond(Matrix const &SK, std::vector<double> const &r, std::vector<double> &w,
                 double const omega=1.0);



/**
 * @brief The Multigrid hierarchy including meshes, vectors and matrices, prolongations is stored.
*/
class Multigrid
{
    public:
       /**
		 * Generates the mesh hierachy with @p nlevel meshes starting from coarse mesh @p cmesg .
         *
         * The refined meshes are generated by edge bisection.
         * All memory is allocated but stiffness matrices are yet not calculated
		 *
		 * @param[in] cmesh	  initial coarse mesh
		 * @param[in] nlevel  number of meshes in hierarchy, including the initial coarse mesh
		 *
		*/
       Multigrid(Mesh const& cmesh, int nlevel);

       Multigrid(Multigrid const&)            = delete;
       Multigrid& operator=(Multigrid const&) = delete;

       ~Multigrid();

       /**
		 * @return  Number of meshes in hierarchy.
		 */
       size_t Nlevels() const
       {return _meshes.size(); }

       /**
		 * @return  Number of Unknowns.
		 */
       int Ndofs() const
       {return _meshes[Nlevels()-1].Nnodes(); }

       /**
		 * @return  Meshes number @p lev .
		*/
       Mesh const& GetMesh(int lev) const
       { return _meshes[lev]; }

       /**
		 * @return  Solution vector at level @p lev .
		 */
       std::vector<double> const&  GetSolution(int lev) const
       { return _u.at(lev); }

       /**
		 * Calculates PDE matrices for all levels.
		 */
       void DefineOperators();

       /**
		 * Calculates PDE matrix for level @p lev.
         *
         * @param[in] lev  level in hierachy
		 */
       void DefineOperator(size_t lev);

       /**
		 * Solves the system of equations at level @p lev via Jacobi iterations
         *
         * @param[in] lev  level in hierachy
		 */
       void JacobiSolve(size_t lev);

      /**
		 * Peformes one multigrid step at level @p lev .
         *
         * @param[in] lev         level in hierachy
         * @param[in] pre_smooth  number of pre/post-smoothing sweeps
         * @param[in] bzero       start with zero-vector as solution
         * @param[in] nu          defines the multigrid cycle (1/2 = V/W)
		 */
       void MG_Step(size_t lev, int pre_smooth=1, bool const bzero=false, int nu=1);

      /**
		 * Solves the system of equations at finest level via multigrid
         *
         * @param[in] pre_smooth  number of pre/post-smoothing sweeps
         * @param[in] eps         stopping criteria (relative error)
         * @param[in] nu          defines the multigrid cycle (1/2 = V/W)
		 */
       void MG_Solve(int pre_smooth=1, double eps=1e-6, int nu=1);


    private:
       gMesh_Hierarchy _meshes;                //!< mesh hierarchy from coarse (level 0) to fine.
       std::vector<FEM_Matrix>           _SK;  //!< Sparse matrix on each level.
       std::vector<std::vector<double>>  _u;   //!< Solution vector on each level.
       std::vector<std::vector<double>>  _f;   //!< Right hand side vector on each level.
       std::vector<std::vector<double>>  _d;   //!< Defect vector on each level.
       std::vector<std::vector<double>>  _w;   //!< Correction vector on each level.
       std::vector<BisectIntDirichlet> _Pc2f;  //!< Interpolation to level from next coarser level.

};




#endif