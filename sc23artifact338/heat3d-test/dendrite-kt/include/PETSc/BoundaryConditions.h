//
// Created by maksbh on 11/28/18.
//

#ifndef DENDRITE2_0_BOUNDARYCONDITIONS_H
#define DENDRITE2_0_BOUNDARYCONDITIONS_H


#include <petscmat.h>
#include <oda.h>
#include <DataTypes.h>
#include <talyfem/grid/zeroptv.h>
#include <OctToPhysical.h>

namespace PETSc{
    typedef std::pair<int, PetscScalar> Dirichlet;

    struct Boundary {
        std::vector<Dirichlet> dirichlets;  // dof to value

        inline void addDirichlet(int dof, PetscScalar value) {
            for (unsigned int i = 0; i < dirichlets.size(); i++) {
                if (dirichlets[i].first == dof) {
                    if (dirichlets[i].second != value)
                        throw std::runtime_error("Cannot set same dirichlet condition multiple times to different values");

                    return;  // dont push back; value already set
                }
            }

            dirichlets.emplace_back(std::make_pair(dof, value));
        }

        inline bool empty() const { return dirichlets.empty(); }
    };

    class BoundaryConditions {
        DENDRITE_UINT ndof_;
        OctToPhysical octToPhys_;


    protected:
        std::vector<PetscInt> m_boundaryRows;
        std::vector<PetscScalar> m_boundaryValues;
    #ifdef IBM
        const std::vector<PetscInt> * ibmDirichletNodes_;
    #endif

    public:

        void clear();

        void init(OctToPhysical & octToPhysical);

        void addByNodalFunction(ot::DA<DIM>* da, const int ndof,
                            const std::function<Boundary(TALYFEMLIB::ZEROPTV, std::size_t nodeID)> &f);

//        void addDirichlet(PetscInt row, PetscScalar val);

        PetscErrorCode applyMatBC(ot::DA<DIM> *da, Mat mat);

//        PetscErrorCode applyMatIBMBC(ot::DA<DIM> *da, Mat mat, const std::vector<PetscInt> * dirichletNodeIDs);

        // sets the values in rhs that have dirichlet boundary conditions to their values
        // note: do not call this on the SNES residual - the SNES residual is actually the delta, which should probably be 0
        //       call this on the solution, then use applyResidualBC on the residual

        PetscErrorCode applyVecBC(ot::DA<DIM> *da, Vec rhs);
#ifdef IBM
        inline void assignIBMDirichletNodes(const std::vector<PetscInt> * dirichletNodes){
          ibmDirichletNodes_ = dirichletNodes;
        }
#endif

//        PetscErrorCode applyVecIBMBC(ot::DA<DIM>* da, Vec rhs, const std::vector<PetscInt> * dirichletNodeIDs);
        // zeros the rows in residual that have dirichlet boundary conditions
        // (i.e. sets the delta to zero)

        PetscErrorCode applyResidualBC(ot::DA<DIM> *da, Vec residual);

        // sets out[i] = in[i] for all rows i that have dirichlet boundary conditions (i.e. sets dirichlet rows to diagonal 1)
        PetscErrorCode applyMatrixFreeBC(ot::DA<DIM>* da, Vec in, Vec out);

        inline std::vector<PetscInt>& rows() {
            return m_boundaryRows;
        }
        inline std::vector<PetscScalar>& values() {
            return m_boundaryValues;
        }

        inline void setDof(DENDRITE_UINT ndof){
            ndof_ = ndof;
        }
#ifdef IBM

//        PetscErrorCode fillWithDirichletID_IBM(DA* da,Mat mat);
//        PetscErrorCode applyIBMBoundaryCondition(DA* da, Mat mat);
//        PetscErrorCode applyIBMBoundaryCondition(DA* da, Mat mat, const std::vector<PetscInt> * dirichletNodeIDs);

#endif
    };
}
#endif //DENDRITE2_0_BOUNDARYCONDITIONS_H
