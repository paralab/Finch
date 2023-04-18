//
// Created by milinda on 12/25/16.
//

/**
 * @author Milinda Fernando
 * @author Hari Sundar
 * @breif Contains data structures to store the reference element information.
 *
 * @refference: Based of HOMG code written in matlab.
 * */
#ifndef SFCSORTBENCH_REFERENCEELEMENT_H
#define SFCSORTBENCH_REFERENCEELEMENT_H

#include "basis.h"
#include "tensor.h"
#include "mathUtils.h"
#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include "interpMatrices.h"
template <typename T>
void printArray_1D(const T *a, int length)
{
    for (int i = 0; i < length; i++)
    {
        std::cout << a[i] << " ";
    }
    std::cout << std::endl;
}

template <typename T>
void printArray_2D(const T *a, int length1, int length2)
{
    for (int i = 0; i < length1; i++)
    {
        for (int j = 0; j < length2; j++)
        {
            std::cout << a[i * length2 + j] << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

class RefElement
{

  private:
    /** Dimension */
    int m_uiDimension;
    /** Polynomial Order */
    int m_uiOrder;
    /** Number of 3D interpolation points on the element */
    int m_uiNp;
    /// /** Number of 2D face interpolation points */
    /// int m_uiNfp;
    /** Number of 1D interpolation points */
    int m_uiNrp;

    bool m_isValid;

    /** reference element volume */
    unsigned int m_uiVol;

    /** 1D reference coordinates of the interpolation nodes (uniform nodal points) */
    std::vector<double> u;

    /** 1D reference coordinates of the interpolation nodes (gll points)*/
    std::vector<double> r;

    /** 1D regular points corresponding to child 0 of u*/
    std::vector<double> u_0;

    /** 1D regular points corresponding to child 1 of u*/
    std::vector<double> u_1;

    /** 1D Gauss points (used for quadrature)*/
    std::vector<double> g;

    /** 1D weights for gauss quadrature */
    std::vector<double> w;

    /** 1D weights for gll quadrature*/
    std::vector<double> wgll;

    /** 1D interpolation matrix for child 0 */
    std::vector<double> ip_1D_0;

    /** 1D interpolation matrix for child 1*/
    std::vector<double> ip_1D_1;

    /** 1D interpolation matrix for child 0 (transpose) */
    std::vector<double> ipT_1D_0;

    /** 1D interpolation matrix for child 1 (transpose)*/
    std::vector<double> ipT_1D_1;

    /**Vandermonde matrix for interpolation points r.   */
    std::vector<double> Vr;

    /**Vandermonde matrix for interpolation points u.   */
    std::vector<double> Vu;

    /**Vandermonde matrix for polynomial at gauss points */
    std::vector<double> Vg;

    /**gradient of the vandermonde for polynomial eval at points u*/
    std::vector<double> gradVu;

    /**gradient of the vandermonde for polynomial eval at points r*/
    std::vector<double> gradVr;

    /**gradient of the vandermonde for polynomial eval at points g*/
    std::vector<double> gradVg;

    /**derivative of the pol. eval at points r. */
    std::vector<double> Dr;

    /** derivative of the pol. eval at the gauss points. */
    std::vector<double> Dg;

    /** derivative of the pol. eval at the gauss points. (transpose) */
    std::vector<double> DgT;

    /** 1D quadrature matrix*/
    std::vector<double> quad_1D;

    /** 1D quadrature matrix transpose*/
    std::vector<double> quadT_1D;


    // Could avoid storing the entrywise squares if simplify KroneckerProduct
    // to allow squaring on the fly.

    /** Entrywise square of quadT_1D, used in computing elemental diagonals. */
    std::vector<double> quadT_1D_hadm2;

    /** Entrywise square of DgT, used in computing elemental diagonals. */
    std::vector<double> DgT_hadm2;


    /**Vandermonde matrix for interpolation points of child 0   */
    std::vector<double> Vu_0;

    /**Vandermonde matrix for interpolation points of child 1   */
    std::vector<double> Vu_1;

    /**intermediate vec 1 needed during interploation */
    std::vector<double> im_vec1;

    /**intermediate vec 1 needed during interploation */
    std::vector<double> im_vec2;


    template <unsigned int dim>   // TODO make this a class template parameter
    inline void getIpPtrs(const double * ipAxis[], unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      const double * const ip[2] = {&(*ip_1D_0.cbegin()), &(*ip_1D_1.cbegin())};
      #pragma unroll(dim)
      for (int d = 0; d < dim; d++)
        ipAxis[d] = ip[(bool) (childNum & (1u<<d))];
    }

    template <unsigned int dim>   // TODO make this a class template parameter
    inline void getIpTPtrs(const double * ipTAxis[], unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      const double * const ipT[2] = {&(*ipT_1D_0.cbegin()), &(*ipT_1D_1.cbegin())};
      #pragma unroll(dim)
      for (int d = 0; d < dim; d++)
        ipTAxis[d] = ipT[(bool) (childNum & (1u<<d))];
    }

    template <typename da, unsigned int dim>
    inline void getDoubleBufferPipeline(const da * fromPtrs[], da * toPtrs[], const da * in, da * out) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      toPtrs[0] = getImVec1();
      fromPtrs[0] = getImVec2();
      for (int d = 1; d < dim; d++)
      {
        fromPtrs[d] = const_cast<const da *>(toPtrs[d-1]);
        toPtrs[d] = const_cast<da *>(fromPtrs[d-1]);
      }
      fromPtrs[0] = in;
      toPtrs[dim-1] = out;
    }



  public:
    RefElement();
    RefElement(unsigned int dim, unsigned int order);
    ~RefElement();
    RefElement(const RefElement &other) = delete;
    RefElement & operator= (const RefElement &other) = delete;

    RefElement(RefElement &&other);
    RefElement & operator= (RefElement &&other);

    // some getter methods to access required data.
    inline int getOrder() const { return m_uiOrder; }
    inline int getDim() const { return m_uiDimension; }
    inline int get1DNumInterpolationPoints() const { return m_uiNrp; }

    inline const double *getIMChild0() const { if (!m_isValid) throw "Invalid!"; return &(*(ip_1D_0.begin())); }
    inline const double *getIMChild1() const { if (!m_isValid) throw "Invalid!"; return &(*(ip_1D_1.begin())); }

    inline const double *getQ1d() const { if (!m_isValid) throw "Invalid!"; return &(*(quad_1D.begin())); }
    inline const double *getQT1d() const { if (!m_isValid) throw "Invalid!"; return &(*(quadT_1D.begin())); }
    inline const double *getDg1d() const { if (!m_isValid) throw "Invalid!"; return &(*(Dg.begin())); }
    inline const double *getDgT1d() const { if (!m_isValid) throw "Invalid!"; return &(*(DgT.begin())); }
    inline const double *getDr1d() const { if (!m_isValid) throw "Invalid!"; return &(*(Dr.begin())); }

    inline const double *getQT1d_hadm2()  const { if (!m_isValid) throw "Invalid!"; return &(*quadT_1D_hadm2.begin()); }
    inline const double *getDgT1d_hadm2() const { if (!m_isValid) throw "Invalid!"; return &(*DgT_hadm2.begin()); }

    inline double *getImVec1() const { if (!m_isValid) throw "Invalid!"; return (double *) &(*(im_vec1.begin())); }
    inline double *getImVec2() const { if (!m_isValid) throw "Invalid!"; return (double *) &(*(im_vec2.begin())); }

    inline const double *getWgq() const { if (!m_isValid) throw "Invalid!"; return &(*(w.begin())); }
    inline const double *getWgll() const { if (!m_isValid) throw "Invalid!"; return &(*(wgll.begin())); }

    inline double getElementSz() const { if (!m_isValid) throw "Invalid!"; return (u.back() - u.front()); }

     /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: interpolated values.
     *
     * @brief This is computed in way that 4d coordinates changes in the order of t,z, y, x
     * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
     *
     * @note Computations are done in internal buffers. It is safe to re-use out == in.
     */
    template <unsigned int dim>
    inline void IKD_Parent2Child(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      assert((childNum < (1u<<dim)));

      // Double buffering.
      const double * imFrom[dim];
      double * imTo[dim];
      getDoubleBufferPipeline<double, dim>(imFrom, imTo, in, out);
      if (dim == 1 && in == out)
        imTo[0] = getImVec1();   // Protect 'in'.

      // Line up 1D operators for each axis, based on childNum.
      const double * ipAxis[dim];
      getIpPtrs<dim>(ipAxis, childNum);

      // Hard-code for now. TODO remove this, it's just here as a backup.
      /// IterateTensorBindMatrix<dim, 0>::template iterate_bind_matrix<double>(m_uiNrp, ipAxis[0], imFrom[0], imTo[0]);
      /// IterateTensorBindMatrix<dim, 1>::template iterate_bind_matrix<double>(m_uiNrp, ipAxis[1], imFrom[1], imTo[1]);
      /// IterateTensorBindMatrix<dim, 2>::template iterate_bind_matrix<double>(m_uiNrp, ipAxis[2], imFrom[2], imTo[2]);
      /// IterateTensorBindMatrix<dim, 3>::template iterate_bind_matrix<double>(m_uiNrp, ipAxis[3], imFrom[3], imTo[3]);

      KroneckerProduct<dim, double, true>(m_uiNrp, ipAxis, imFrom, imTo, ndofs);

      if (dim == 1 && in == out)   // Protected 'in'.
        memcpy(out, imTo[dim-1], sizeof(double)*ndofs*intPow(m_uiNrp, dim));
    }

    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the contribution needed to be computed.
     * @param[out] out: child to parent contribution values. (used in FEM integral ivaluation)
     *
     * @brief This is computed in way that 3d coordinates changes in the order of z, y, x
     * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
     *
     */
    template <unsigned int dim>
    inline void IKD_Child2Parent(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      assert((childNum < (1u<<dim)));

      // Double buffering.
      const double * imFrom[dim];
      double * imTo[dim];
      getDoubleBufferPipeline<double, dim>(imFrom, imTo, in, out);
      if (dim == 1 && in == out)
        imTo[0] = getImVec1();   // Protect 'in'.

      // Line up 1D operators for each axis, based on childNum.
      const double * ipTAxis[dim];
      getIpTPtrs<dim>(ipTAxis, childNum);

      KroneckerProduct<dim, double, true>(m_uiNrp, ipTAxis, imFrom, imTo, ndofs);

      if (dim == 1 && in == out)   // Protected 'in'.
        memcpy(out, imTo[dim-1], sizeof(double)*ndofs*intPow(m_uiNrp, dim));
    }


    // ------------------------------------------- //



    inline void I4D_Parent2Child(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      IKD_Parent2Child<4>(in, out, ndofs, childNum);
    }

    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: interpolated values.
     *
     * @brief This is computed in way that 3d coordinates changes in the order of z, y, x
     * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
     *
     */
    inline void I3D_Parent2Child(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Parent2Child<3>(in, out, ndofs, childNum);

        /// double *im1 = (double *)&(*(im_vec1.begin()));
        /// double *im2 = (double *)&(*(im_vec2.begin()));

        /// switch (childNum)
        /// {
        /// case 0:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 1:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im2, out); // along z

        ///     break;
        /// case 2:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 3:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 4:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 5:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 6:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 7:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ip_1D_1.begin())), im2, out); // along z
        ///     break;
        /// default:
        ///     std::cout << "[refel][error]: invalid child number specified for 3D interpolation." << std::endl;
        ///     break;
        /// }
    }



    inline void I4D_Child2Parent(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

      IKD_Child2Parent<4>(in, out, ndofs, childNum);
    }



    /**
    * @param[in] in: input function values.
    * @param[in] childNum: Morton ID of the child number where the contribution needed to be computed.
    * @param[out] out: child to parent contribution values. (used in FEM integral ivaluation)
    *
    * @brief This is computed in way that 3d coordinates changes in the order of z, y, x
    * Which means first we need to fill all the z values in plane(x=0,y=0) then all the z values in plane (x=0,y=0+h) and so forth.
    *
    */

    inline void I3D_Child2Parent(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Child2Parent<3>(in, out, ndofs, childNum);

        /// double *im1 = (double *)&(*(im_vec1.begin()));
        /// double *im2 = (double *)&(*(im_vec2.begin()));

        /// switch (childNum)
        /// {
        /// case 0:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 1:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im2, out); // along z

        ///     break;
        /// case 2:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 3:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im2, out); // along z
        ///     break;
        /// case 4:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 5:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 6:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im2, out); // along z
        ///     break;
        /// case 7:
        ///     DENDRO_TENSOR_IIAX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_IAIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im1, im2); // along y
        ///     DENDRO_TENSOR_AIIX_APPLY_ELEM(m_uiNrp, &(*(ipT_1D_1.begin())), im2, out); // along z
        ///     break;
        /// default:
        ///     std::cout << "[refel][error]: invalid child number specified for 3D interpolation." << std::endl;
        ///     break;
        /// }

#ifdef FEM_ACCUMILATE_ONES_TEST
        for (unsigned int node = 0; node < (m_uiNrp * m_uiNrp * m_uiNrp); node++)
            out[node] = 1.0;
#endif
    }

    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: interpolated values.
     * */

    inline void I2D_Parent2Child(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Parent2Child<2>(in, out, ndofs, childNum);

        /// double *im1 = (double *)&(*(im_vec1.begin()));
        /// double *im2 = (double *)&(*(im_vec2.begin()));

        /// switch (childNum)
        /// {

        /// case 0:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_0.begin())), im1, out); // along y (in 3d z)
        ///     break;
        /// case 1:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_0.begin())), im1, out); // along y (in 3d z)
        ///     break;

        /// case 2:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_1.begin())), im1, out); // along y (in 3d z)
        ///     break;

        /// case 3:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ip_1D_1.begin())), im1, out); // along y (in 3d z)
        ///     break;
        /// default:
        ///     std::cout << "[refel][error]: invalid child number specified for 2D  interpolation." << std::endl;
        ///     break;
        /// }
    }

    /**
     * @param[in] in: input function values.
     * @param[in] childNum: Morton ID of the child number where the interpolation is needed.
     * @param[out] out: child to parent contribution values. (used in FEM integral ivaluation)
     * */

    inline void I2D_Child2Parent(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Child2Parent<2>(in, out, ndofs, childNum);

        /// double *im1 = (double *)&(*(im_vec1.begin()));
        /// double *im2 = (double *)&(*(im_vec2.begin()));

        /// switch (childNum)
        /// {

        /// case 0:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_0.begin())), im1, out); // along y (in 3d z)
        ///     break;
        /// case 1:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_0.begin())), im1, out); // along y (in 3d z)
        ///     break;

        /// case 2:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_0.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_1.begin())), im1, out); // along y (in 3d z)
        ///     break;

        /// case 3:
        ///     DENDRO_TENSOR_IAX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_1.begin())), in, im1);  // along x
        ///     DENDRO_TENSOR_AIX_APPLY_ELEM_2D(m_uiNrp, &(*(ipT_1D_1.begin())), im1, out); // along y (in 3d z)
        ///     break;
        /// default:
        ///     std::cout << "[refel][error]: invalid child number specified for 2D  interpolation." << std::endl;
        ///     break;
        /// }

#ifdef FEM_ACCUMILATE_ONES_TEST
        for (unsigned int node = 0; node < (m_uiNrp * m_uiNrp); node++)
            out[node] = 1.0;
#endif
    }

    /**
     * @param [in] in input function values
     * @param [in] childNum Morton ID of the child number
     * @param [out] interpolated values from parent to child.
     * */

    inline void I1D_Parent2Child(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Parent2Child<1>(in, out, ndofs, childNum);

        /// switch (childNum)
        /// {
        /// case 0:
        ///     for (unsigned int i = 0; i < m_uiNrp; i++)
        ///     {
        ///         out[i] = 0.0;
        ///         for (unsigned int j = 0; j < m_uiNrp; j++)
        ///         {
        ///             out[i] += ip_1D_0[j * m_uiNrp + i] * in[j];
        ///         }
        ///     }
        ///     break;
        /// case 1:
        ///     for (unsigned int i = 0; i < m_uiNrp; i++)
        ///     {
        ///         out[i] = 0.0;
        ///         for (unsigned int j = 0; j < m_uiNrp; j++)
        ///         {
        ///             out[i] += ip_1D_1[j * m_uiNrp + i] * in[j];
        ///         }
        ///     }
        ///     break;

        /// default:
        ///     std::cout << "[refel][error]: Invalid child number specified for 1D interpolation. " << std::endl;
        ///     break;
        /// }
    }

    /**
     * @param [in] in input function values
     * @param [in] childNum Morton ID of the child number
     * @param [out] child to parent contribution values. (used in FEM integral ivaluation)
     * */

    inline void I1D_Child2Parent(const double *in, double *out, unsigned int ndofs, unsigned int childNum) const
    {
      if (!m_isValid)
        throw "RefElement constructed or copied incorrectly!";

        IKD_Child2Parent<1>(in, out, ndofs, childNum);

        /// switch (childNum)
        /// {
        /// case 0:
        ///     for (unsigned int i = 0; i < m_uiNrp; i++)
        ///     {
        ///         out[i] = 0.0;
        ///         for (unsigned int j = 0; j < m_uiNrp; j++)
        ///         {
        ///             out[i] += ipT_1D_0[j * m_uiNrp + i] * in[j];
        ///         }
        ///     }
        ///     break;
        /// case 1:
        ///     for (unsigned int i = 0; i < m_uiNrp; i++)
        ///     {
        ///         out[i] = 0.0;
        ///         for (unsigned int j = 0; j < m_uiNrp; j++)
        ///         {
        ///             out[i] += ipT_1D_1[j * m_uiNrp + i] * in[j];
        ///         }
        ///     }
        ///     break;

        /// default:
        ///     std::cout << "[refel][error]: Invalid child number specified for 1D interpolation. " << std::endl;
        ///     break;
        /// }

#ifdef FEM_ACCUMILATE_ONES_TEST
        for (unsigned int node = 0; node < (m_uiNrp); node++)
            out[node] = 1.0;
#endif
    }

    void generateHeaderFile(char *fName);
};


namespace detail
{
  template <typename NodeT>
  struct VectorCast
  {
    template <typename InputT>
    static
    std::vector<NodeT> create(const InputT *begin, size_t size)
    {
      std::vector<NodeT> result(begin, begin + size);
      return result;
    }
  };

  template <typename NodeT>
  struct VectorTest
  {
    template <typename InputT>
    static
    std::vector<NodeT> create(const InputT *begin, size_t size)
    {
      std::vector<NodeT> result;
      for (size_t ii = 0; ii < size; ++ii)
        result.push_back(fabs(begin[ii]) > 1e-5);
      return result;
    }
  };

  template <typename NodeT>
  struct CastFloatsTestInts : public VectorTest<NodeT> { };

  template <>
  struct CastFloatsTestInts<float> : public VectorCast<float> { };
  template <>
  struct CastFloatsTestInts<double> : public VectorCast<double> { };
  template <>
  struct CastFloatsTestInts<long double> : public VectorCast<long double> { };
}

template <unsigned int dim, typename NodeT>
struct InterpMatrices
{
  std::vector<NodeT> m_ip_1D[2];   // Parent to child interpolation.
  std::vector<NodeT> m_ipT_1D[2];  // Child to parent.
  mutable std::vector<NodeT> m_imBufs[2];    // Intermediate buffers.

  unsigned int m_eleOrder;
  unsigned int m_npe;

  InterpMatrices() : InterpMatrices(1) {}

  InterpMatrices(unsigned int eleOrder)
  {
    m_eleOrder = eleOrder;
    m_npe = intPow(eleOrder+1, dim);

    const unsigned int ipMatSz = (eleOrder+1)*(eleOrder+1);
    RefElement tempRefEl(dim, eleOrder);
    m_ip_1D[0] = detail::CastFloatsTestInts<NodeT>::create(tempRefEl.getIMChild0(), ipMatSz);
    m_ip_1D[1] = detail::CastFloatsTestInts<NodeT>::create(tempRefEl.getIMChild1(), ipMatSz);
    m_ipT_1D[0].resize(ipMatSz);
    m_ipT_1D[1].resize(ipMatSz);
    for (int ii = 0; ii < eleOrder+1; ii++)     // Transpose
      for (int jj = 0; jj < eleOrder+1; jj++)
      {
        m_ipT_1D[0][ii * (eleOrder+1) + jj] = m_ip_1D[0][jj * (eleOrder+1) + ii];
        m_ipT_1D[1][ii * (eleOrder+1) + jj] = m_ip_1D[1][jj * (eleOrder+1) + ii];
      }
  }

  static constexpr bool P2C = false;
  static constexpr bool C2P = true;

  template <bool useTranspose>
  inline void IKD_ParentChildInterpolation(const NodeT *in,
                                           NodeT *out,
                                           unsigned int ndofs,
                                           unsigned int childNum_m) const
  {
    assert(dim > 1 || in != out);  // Not enough intermediate buffers.

    const std::vector<NodeT> (&ip_1D)[2] = (useTranspose ? m_ipT_1D : m_ip_1D);

    // Line up 1D operators for each axis, based on childNum.
    const NodeT *ipAxis[dim];
    for (int d = 0; d < dim; d++)
      ipAxis[d] = ip_1D[bool(childNum_m & (1u << d))].data();

    // Double buffering of parent node coordinates during interpolation.
    const NodeT * imFrom[dim];
    NodeT * imTo[dim];
    m_imBufs[0].resize(ndofs * m_npe);
    m_imBufs[1].resize(ndofs * m_npe);
    for (int d = 0; d < dim; d++)
    {
      imTo[d] = &(*m_imBufs[d % 2].begin());
      imFrom[d] = &(*m_imBufs[!(d % 2)].begin());
    }
    imFrom[0] = in;     // Overwrite pointer to first source.
    imTo[dim-1] = out;  // Overwrite pointer to final destination.

    // Interpolate all element nodes.
    // (The ones we actually use should have valid values.)
    KroneckerProduct<dim, NodeT, true>(m_eleOrder+1, ipAxis, imFrom, imTo, ndofs);
  }
};

#endif //SFCSORTBENCH_REFERENCEELEMENT_H
