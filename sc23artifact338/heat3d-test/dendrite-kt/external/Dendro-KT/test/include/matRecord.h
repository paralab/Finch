//
// Created by milinda on 12/20/18.
//

/**
 * @brief Matrix record class to be used while setting values in the global matrix.
 * This is based from the old dendro MatRecord class, where now it supports  both petsc and non petsc matrices.
 *
 * @author Milinda Fernando, School of Computing, University of Utah.
 * */

#ifdef BUILD_WITH_PETSC
    #include "petsc.h"
#endif

#include <iostream>
#include <ostream>
#include "dendro.h"

#ifndef DENDRO_KT_MATRECORD_H
#define DENDRO_KT_MATRECORD_H

namespace ot
{

    class MatRecord {

        private:
            /**@brief local index of the node coresponding to the local row of the matrix */
            RankI m_uiRowID;
            /**@brief local index of the node corresponding to the local column of the matrix */
            RankI m_uiColID;
            /**@brief The degree of freedom (0-based) of the node corresponding to  the row of the matrix.*/
            RankI m_uiRowDim;
            /**@brief The degree of freedom (0-based) of the node corresponding to  the column of the matrix. */
            RankI m_uiColDim;

#ifdef BUILD_WITH_PETSC
            /**@brief  matrix value (petsc)*/
            PetscScalar m_uiVal;
#else
            /**@brief  matrix value (dendro)*/
            DendroScalar m_uiVal;
#endif

        public:

            /**
             *@brief The default constructor
             */
            MatRecord() {
                m_uiRowID = m_uiColID = m_uiRowDim = m_uiColDim = static_cast<RankI>(-1);
                m_uiVal = 0;
            }

            /**@brief returns the row ID*/
            inline RankI getRowID() const {return m_uiRowID;}
            /**@brief returns the col ID*/
            inline RankI getColID() const {return m_uiColID;}
            /**@brief return the row dof*/
            inline RankI getRowDim() const {return m_uiRowDim;}
            /**@brief return the col dof*/
            inline RankI getColDim() const {return m_uiColDim;}

            /**@brief sets the rowID value*/
            inline void setRowID(RankI rowID) {m_uiRowID=rowID;}

            /**@brief sets the rowID value*/
            inline void setColID(RankI colID) {m_uiColID=colID;}

            /**@brief sets the rowDim value*/
            inline void setRowDim(RankI rowDim) {m_uiRowDim=rowDim;}

            /**@brief sets the colDim value*/
            inline void setColDim(RankI colDim) {m_uiColDim=colDim;}


#ifdef BUILD_WITH_PETSC

            /**@breif MatRecord constructor
             * @param [in] rowID: row id
             * @param [in] colID: col id
             * @param [in] rowDim: row dim
             * @param [in] colDim: colDim
             * @param [in] value: value of for the entry.
             * */

            MatRecord(RankI rowID, RankI colID, RankI rowDim, RankI colDim, PetscScalar value)
            {
                m_uiRowID=rowID;
                m_uiColID=colID;

                m_uiRowDim=rowDim;
                m_uiColDim=colDim;

                m_uiVal=value;

            }



            /**@brief returns the entry value*/
            inline PetscScalar getMatVal() const {return m_uiVal;}

            /**@brief sets matrix value*/
            inline void setMatValue(PetscScalar value) {m_uiVal=value;}

#else

            /**@breif MatRecord constructor
             * @param [in] rowID: row id
             * @param [in] colID: col id
             * @param [in] rowDim: row dim
             * @param [in] colDim: colDim
             * @param [in] value: value of for the entry.
             * */

            MatRecord(RankI rowID, RankI colID, RankI rowDim, RankI colDim, DendroScalar value)
            {
                m_uiRowID=rowID;
                m_uiColID=colID;

                m_uiRowDim=rowDim;
                m_uiColDim=colDim;

                m_uiVal=value;

            }



            /**@brief returns the entry value*/
            inline DendroScalar getMatVal() const {return m_uiVal;}

            /**@brief sets matrix value*/
            inline void setMatValue(DendroScalar value) {m_uiVal=value;}

#endif

            /**@brief The copy constructor */
            MatRecord(const MatRecord & other) {
                m_uiRowID = other.getRowID();
                m_uiColID = other.getColID();
                m_uiRowDim = other.getRowDim();
                m_uiColDim = other.getColDim();
                m_uiVal = other.getMatVal();
            }


            /** @brief The assignment operator */
            MatRecord & operator = (MatRecord const  & other) {
                if(this == (&other)) {return *this;}
                m_uiRowID = other.getRowID();
                m_uiColID = other.getColID();
                m_uiRowDim = other.getRowDim();
                m_uiColDim = other.getColDim();
                m_uiVal = other.getMatVal();
                return *this;
            }

            /** @brief Overloaded == Operator */
            bool  operator == ( MatRecord const  &other) const {
                return ( (m_uiRowID == other.getRowID()) && (m_uiColID == other.getColID())
                         && (m_uiVal == other.getMatVal()) && (m_uiRowDim == other.getRowDim()) && (m_uiColDim == other.getColDim()) );
            }

            /** @brief Overloaded != Operator */
            bool  operator != (MatRecord const  &other) const {
                return (!((*this) == other));
            }

            /** @brief Overloaded < Operator */
            bool  operator < (MatRecord const  &other) const {
                if(m_uiRowID < other.getRowID()) {
                    return true;
                }else if(m_uiRowID == other.getRowID()) {
                    if(m_uiRowDim < other.getRowDim()) {
                        return true;
                    }else if(m_uiRowDim == other.getRowDim()) {
                        if(m_uiColID < other.getColID()) {
                            return true;
                        }else if (m_uiColID == other.getColID()) {
                            if(m_uiColDim < other.getColDim()) {
                                return true;
                            }else if(m_uiColDim == other.getColDim()) {
                                return (m_uiVal < other.getMatVal());
                            }else {
                                return false;
                            }
                        }else {
                            return false;
                        }
                    }else {
                        return false;
                    }
                }else {
                    return false;
                }
            }

            /** @brief Overloaded > Operator */
            bool  operator > (MatRecord const  &other) const {
                return ( (!((*this) < other)) && ((*this) != other) );
            }

            /** @brief Overloaded <= Operator */
            bool  operator <= (MatRecord const  &other) const {
                return ( ((*this) < other) || ((*this) == other) );
            }

            /** @brief Overloaded >= Operator */
            bool  operator >= (MatRecord const  &other) const {
                return (!((*this) < other)) ;
            }

            friend std::ostream & operator<< (std::ostream & os, MatRecord const & re) 
            {
                return (os << " row : "<<re.getRowID() << " col: "<<re.getColID()<<" rdim: "<<re.getRowDim()<<" cdim: "<<re.getColDim() );
            }

            

    };



    /**
     * MatCompactRows
     *
     * @brief Store entries of the matrix in chunks corresponding to rows of elemental matrices.
     *
     * There are (npe x ndofs) entries per chunk.
     * Multiple instances per entry are allowed. They will get summed later.
     * Storage ordering assumes dofs index faster than nodes.
     * Row and column indices have degrees of freedom built in.
     * That is, colIdx==(colNode*ndofs + dofIdx).
     * Same convention applies to appendChunk().
     */
    class MatCompactRows
    {
      public:
#ifdef BUILD_WITH_PETSC
        using ScalarT = PetscScalar;
        using IndexT = PetscInt;
#else
        using ScalarT = DendroScalar;
        using IndexT = DendroIntL;
#endif

      protected:
        unsigned int m_ndofs;
        unsigned int m_npe;
        std::vector<IndexT> m_rowIdxs;
        std::vector<IndexT> m_colIdxs;  // size==m_rowIdxs.size()*m_npe*m_ndofs.
        std::vector<ScalarT> m_colVals; // size==m_rowIdxs.size()*m_npe*m_ndofs.


      public:
        MatCompactRows() = delete;
        MatCompactRows(unsigned int npe, unsigned int ndofs)
          : m_npe(npe), m_ndofs(ndofs)
        {}

        // Should probably leave copy constructor / copy assignment deleted.
        // Should probably define move constructor / move assignment if needed.

        //
        // Getters.
        //
        unsigned int getNdofs()     const { return m_ndofs; }
        unsigned int getNpe()       const { return m_npe; }
        unsigned int getChunkSize() const { return m_npe * m_ndofs; }
        size_t       getNumRows()   const { return m_rowIdxs.size(); }
        const std::vector<IndexT> & getRowIdxs()  const { return m_rowIdxs; }
        const std::vector<IndexT> & getColIdxs()  const { return m_colIdxs; }
        const std::vector<ScalarT> & getColVals() const { return m_colVals; }

        //
        // Mutators.
        //

        /**
         * @brief Add (offset * m_ndofs) to all row and column indices.
         *        E.g. offsetNodeIdxs( da.getGlobalRankBegin() );
         */
        void offsetNodeIdxs(IndexT offset)
        {
          offset *= m_ndofs;
          for (IndexT & ridx : m_rowIdxs)
            ridx += offset;
          for (IndexT & cidx : m_colIdxs)
            cidx += offset;
        }

        /**
         * @brief Append the next row from an elemental matrix.
         */
        template <typename T>
        void appendChunk(IndexT rowIdx, const IndexT *colIdxs, const T *colVals)
        {
          // Force indices to be IndexT so that -1 is -1 (see petsc MatSetValues).

          m_rowIdxs.push_back(rowIdx);
          for (int j = 0; j < m_npe * m_ndofs; j++)
          {
            m_colIdxs.push_back(colIdxs[j]);
            m_colVals.push_back(colVals[j]);
          }
        }

        void clear()
        {
          m_rowIdxs.clear();
          m_colIdxs.clear();
          m_colVals.clear();
        }
    };


}





#endif //DENDRO_KT_MATRECORD_H
