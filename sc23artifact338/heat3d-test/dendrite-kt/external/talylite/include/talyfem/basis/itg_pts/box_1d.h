/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/
// --- end license text --- //
#pragma once


#include <talyfem/grid/zeroptv.h>
#include <talyfem/data_structures/constexpr_array.h>
#include <talyfem/basis/constants.h>

namespace TALYFEMLIB {

#ifdef CLENSHAW_CURTIS
    // TODO: Update the Clenshaw Curtis point.
/**
 * 1D gauss points for N = 1 (1 point).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<1, 1, 0> {
        static constexpr int n_itg_pts = 1;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[1] = {ZEROPTV(0.0, 0, 0)};
        ///! integration point weights
        static constexpr constexpr_array<double, 1> weights = {{2.0}};
    };

    /**
       * 1D Clenshaw-Curtis points for N = 2 (2 points).
   */
    template<>
    struct BoxItgPts<2, 1, 0> {
        static constexpr int n_itg_pts = 2;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[2] = {
                ZEROPTV(-1.0, 0, 0),  // -(1.0 )
                ZEROPTV(+1.0 , 0, 0)   // +(1.0 )
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 2> weights = {{1, 1}};
    };

    /**
 * 1D gauss points for N = 3 (3 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<3, 1, 0> {
        static constexpr int n_itg_pts = 3;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[3] = {
                ZEROPTV(-1.0, 0, 0),  // -(sqrt(3.0) / sqrt(5.0))
                ZEROPTV(0.0, 0, 0),                         // 0
                ZEROPTV(+1.0, 0, 0)   // +(sqrt(3.0) / sqrt(5.0))
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 3> weights =
                {{
                         0.333333333333333333,
                         1.333333333333333333,
                         0.333333333333333333
                 }};
    };

/**
 * 1D gauss points for N = 4 (4 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<4, 1, 0> {
        static constexpr int n_itg_pts = 4;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[4] = {
                ZEROPTV(-1.0, 0, 0),
                ZEROPTV(-0.5, 0, 0),
                ZEROPTV(+0.5, 0, 0),
                ZEROPTV(+1.0, 0, 0)
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 4> weights =
                {{
                         0.111111111111111111,
                         0.888888888888888888,
                         0.888888888888888888,
                         0.111111111111111111
                 }};
    };

/**
 * 1D gauss points for N = 5 (5 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<5, 1, 0> {
        static constexpr int n_itg_pts = 5;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[5] = {
                ZEROPTV(-1.0, 0, 0),
                ZEROPTV(-0.707106781186547573, 0, 0),
                ZEROPTV(0.0, 0, 0),
                ZEROPTV(+0.707106781186547573, 0, 0),
                ZEROPTV(+1.0, 0, 0)
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 5> weights =
                {{
                         0.0666666666666666,
                         0.5333333333333333,
                         0.8,
                         0.5333333333333333,
                         0.0666666666666666
                 }};
    };

/**
 * 1D gauss points for N = 6 (6 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<6, 1, 0> {
        static constexpr int n_itg_pts = 6;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[6] = {
                ZEROPTV(-1.0, 0, 0),
                ZEROPTV(-0.809016994374947451, 0, 0),
                ZEROPTV(-0.30901699437494734, 0, 0),
                ZEROPTV(+0.30901699437494734, 0, 0),
                ZEROPTV(+0.809016994374947451, 0, 0),
                ZEROPTV(+1.0, 0, 0)
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 6> weights =
                {{
                         0.04,
                         0.360743041200011216,
                         0.599256958799988748,
                         0.599256958799988748,
                         0.360743041200011161,
                         0.04
                 }};
    };

/**
 * 1D gauss points for N = 7 (7 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<7, 1, 0> {
        static constexpr int n_itg_pts = 7;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[7] = {
                ZEROPTV(-1.0, 0, 0),  //
                ZEROPTV(-0.866025403784438708, 0, 0),  //
                ZEROPTV(-0.5, 0, 0),  //
                ZEROPTV(0.0, 0, 0),  //
                ZEROPTV(+0.5, 0, 0),  //
                ZEROPTV(+0.866025403784438708, 0, 0),  //
                ZEROPTV(+1.0, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 7> weights =
                {{
                         0.0285714285714285636,
                         0.253968253968254121,
                         0.457142857142857018,
                         0.520634920634920673,
                         0.457142857142857018,
                         0.253968253968254121,
                         0.0285714285714285636
                 }};
    };

/**
 * 1D gauss points for N = 8 (8 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<8, 1, 0> {
        static constexpr int n_itg_pts = 8;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[8] = {
                ZEROPTV(-1.0, 0, 0),  //
                ZEROPTV(-0.900968867902419035, 0, 0),  //
                ZEROPTV(-0.623489801858733483, 0, 0),  //
                ZEROPTV(-0.222520933956314337, 0, 0),  //
                ZEROPTV(+0.222520933956314337, 0, 0),  //
                ZEROPTV(+0.623489801858733483, 0, 0),  //
                ZEROPTV(+0.900968867902419035, 0, 0),  //
                ZEROPTV(+1.0, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 8> weights =
                {{
                         0.020408163265306211,
                         0.190141007218208369,
                         0.352242423718159114,
                         0.437208405798326483,
                         0.437208405798326483,
                         0.352242423718159114,
                         0.190141007218208369,
                         0.020408163265306211
                 }};
    };

/**
 * 1D gauss points for N = 9 (9 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<9, 1, 0> {
        static constexpr int n_itg_pts = 9;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[9] = {
                ZEROPTV(-1.0, 0, 0),  //
                ZEROPTV(-0.923879532511286738, 0, 0),  //
                ZEROPTV(-0.707106781186547573, 0, 0),  //
                ZEROPTV(-0.382683432365089837, 0, 0),  //
                ZEROPTV(0.0, 0, 0),  //
                ZEROPTV(+0.382683432365089837, 0, 0),  //
                ZEROPTV(+0.707106781186547573, 0, 0),  //
                ZEROPTV(+0.923879532511286738, 0, 0),  //
                ZEROPTV(+1.0, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 9> weights =
                {{
                         0.0158730158730158791,
                         0.146218649216018209,
                         0.279365079365079372,
                         0.361717858720489727,
                         0.39365079365079364,
                         0.361717858720489727,
                         0.279365079365079372,
                         0.146218649216018209,
                         0.0158730158730158791
                 }};

    };
    template<>
    struct BoxItgPts<10, 1, 0> {
        static constexpr int n_itg_pts = 10;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[10] = {
                ZEROPTV(-1.0, 0, 0),  //
                ZEROPTV(-0.939692620785908428, 0, 0),  //
                ZEROPTV(-0.766044443118978013, 0, 0),  //
                ZEROPTV(-0.5, 0, 0),  //
                ZEROPTV(-0.173648177666930414, 0, 0),  //
                ZEROPTV(+0.173648177666930414, 0, 0),  //
                ZEROPTV(+0.5, 0, 0),  //
                ZEROPTV(+0.766044443118978013, 0, 0),  //
                ZEROPTV(+0.939692620785908428, 0, 0),
                ZEROPTV(+1.0, 0, 0)//
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 10> weights =
                {{
                         0.0123456790123456263,
                         0.116567456572037209,
                         0.225284323338104553,
                         0.301940035273368546,
                         0.34386250580414407,
                         0.343862505804144458,
                         0.301940035273368546,
                         0.225284323338104553,
                         0.116567456572037209,
                         0.0123456790123456263

             }};
        };
#else
// Gauss point source (Aug 2016):
// https://en.wikipedia.org/wiki/Gaussian_quadrature#Gauss.E2.80.93Legendre_quadrature

// sqrt() cannot be evaluated at compile-time, because it is not constexpr.
// So, these constants are calculated by hand elsewhere.
// The analytic expressions are included in comments where appropriate.

/**
 * 1D gauss points for N = 1 (1 point).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<1, 1, 0> {
        static constexpr int n_itg_pts = 1;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[1] = {ZEROPTV(0.0, 0, 0)};
        ///! integration point weights
        static constexpr constexpr_array<double, 1> weights = {{2.0}};
    };

/**
 * 1D gauss points for N = 2 (2 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<2, 1, 0> {
        static constexpr int n_itg_pts = 2;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[2] = {
                ZEROPTV(-1.0 / 1.7320508075688772935, 0, 0),  // -(1.0 / sqrt(3.0))
                ZEROPTV(+1.0 / 1.7320508075688772935, 0, 0)   // +(1.0 / sqrt(3.0))
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 2> weights = {{1, 1}};
    };

/**
 * 1D gauss points for N = 3 (3 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<3, 1, 0> {
        static constexpr int n_itg_pts = 3;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[3] = {
                ZEROPTV(-0.77459666924148337703585, 0, 0),  // -(sqrt(3.0) / sqrt(5.0))
                ZEROPTV(0.0, 0, 0),                         // 0
                ZEROPTV(+0.77459666924148337703585, 0, 0)   // +(sqrt(3.0) / sqrt(5.0))
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 3> weights =
                {{
                         5.0 / 9.0,
                         8.0 / 9.0,
                         5.0 / 9.0
                 }};
    };

/**
 * 1D gauss points for N = 4 (4 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<4, 1, 0> {
        static constexpr int n_itg_pts = 4;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[4] = {
                ZEROPTV(-0.86113631159405257522394, 0, 0),  // -(sqrt((3/7) + (2/7) * sqrt(6/5)))
                ZEROPTV(-0.33998104358485626480266, 0, 0),  // -(sqrt((3/7) - (2/7) * sqrt(6/5)))
                ZEROPTV(+0.33998104358485626480266, 0, 0),  // +(sqrt((3/7) - (2/7) * sqrt(6/5)))
                ZEROPTV(+0.86113631159405257522394, 0, 0)   // +(sqrt((3/7) + (2/7) * sqrt(6/5)))
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 4> weights =
                {{
                         0.347854845137453857373063,  // (18 - sqrt(30)) / 36
                         0.652145154862546142626936,  // (18 + sqrt(30)) / 36
                         0.652145154862546142626936,  // (18 + sqrt(30)) / 36
                         0.347854845137453857373063   // (18 - sqrt(30)) / 36
                 }};
    };

/**
 * 1D gauss points for N = 5 (5 points).
 * Source: https://en.wikipedia.org/wiki/Gaussian_quadrature
 */
    template<>
    struct BoxItgPts<5, 1, 0> {
        static constexpr int n_itg_pts = 5;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[5] = {
                ZEROPTV(-0.90617984593866399279762, 0, 0),  // -((1/3) * sqrt(5 + 2 * sqrt(10/7)))
                ZEROPTV(-0.53846931010568309103631, 0, 0),  // -((1/3) * sqrt(5 - 2 * sqrt(10/7)))
                ZEROPTV(0.0, 0, 0),                         // 0
                ZEROPTV(+0.53846931010568309103631, 0, 0),  // +((1/3) * sqrt(5 - 2 * sqrt(10/7)))
                ZEROPTV(+0.90617984593866399279762, 0, 0)   // +((1/3) * sqrt(5 + 2 * sqrt(10/7)))
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 5> weights =
                {{
                         0.236926885056189087514264040,  // (322 - 13 * sqrt(70)) / 900
                         0.478628670499366468041291514,  // (322 + 13 * sqrt(70)) / 900
                         128.0 / 225.0,                  // 128.0 / 225.0
                         0.478628670499366468041291514,  // (322 + 13 * sqrt(70)) / 900
                         0.236926885056189087514264040   // (322 - 13 * sqrt(70)) / 900
                 }};
    };

    /**
     * @author maksbh
     * These are added by Robert and copied to TalyLite.
     */

/**
 * 1D gauss points for N = 6 (6 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<6, 1, 0> {
        static constexpr int n_itg_pts = 6;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[6] = {
                ZEROPTV(-0.9324695142031520278123, 0, 0),  //
                ZEROPTV(-0.661209386466264513661, 0, 0),  //
                ZEROPTV(-0.2386191860831969086305, 0, 0),  //
                ZEROPTV(+0.2386191860831969086305, 0, 0),  //
                ZEROPTV(+0.661209386466264513661, 0, 0),  //
                ZEROPTV(+0.9324695142031520278123, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 6> weights =
                {{
                         0.1713244923791703450403,
                         0.3607615730481386075698,
                         0.4679139345726910473899,
                         0.4679139345726910473899,
                         0.3607615730481386075698,
                         0.1713244923791703450403
                 }};
    };

/**
 * 1D gauss points for N = 7 (7 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<7, 1, 0> {
        static constexpr int n_itg_pts = 7;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[7] = {
                ZEROPTV(-0.9491079123427585245262, 0, 0),  //
                ZEROPTV(-0.7415311855993944398639, 0, 0),  //
                ZEROPTV(-0.4058451513773971669066, 0, 0),  //
                ZEROPTV(0.0, 0, 0),  //
                ZEROPTV(+0.4058451513773971669066, 0, 0),  //
                ZEROPTV(+0.7415311855993944398639, 0, 0),  //
                ZEROPTV(+0.9491079123427585245262, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 7> weights =
                {{
                         0.1294849661688696932706,
                         0.2797053914892766679015,
                         0.38183005050511894495,
                         0.417959183673469387755,
                         0.38183005050511894495,
                         0.2797053914892766679015,
                         0.1294849661688696932706
                 }};
    };

/**
 * 1D gauss points for N = 8 (8 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<8, 1, 0> {
        static constexpr int n_itg_pts = 8;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[8] = {
                ZEROPTV(-0.9602898564975362316836, 0, 0),  //
                ZEROPTV(-0.7966664774136267395916, 0, 0),  //
                ZEROPTV(-0.5255324099163289858177, 0, 0),  //
                ZEROPTV(-0.1834346424956498049395, 0, 0),  //
                ZEROPTV(+0.1834346424956498049395, 0, 0),  //
                ZEROPTV(+0.5255324099163289858177, 0, 0),  //
                ZEROPTV(+0.7966664774136267395916, 0, 0),  //
                ZEROPTV(+0.9602898564975362316836, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 8> weights =
                {{
                         0.1012285362903762591525,
                         0.2223810344533744705444,
                         0.313706645877887287338,
                         0.3626837833783619829652,
                         0.3626837833783619829652,
                         0.313706645877887287338,
                         0.2223810344533744705444,
                         0.1012285362903762591525
                 }};
    };

/**
 * 1D gauss points for N = 9 (9 points).
 * Source: https://keisan.casio.com/exec/system/1329114617
 */
    template<>
    struct BoxItgPts<9, 1, 0> {
        static constexpr int n_itg_pts = 9;  ///< number of integration points
        ///! integration points
        static constexpr ZEROPTV itg_pts[9] = {
                ZEROPTV(-0.9681602395076260898356, 0, 0),  //
                ZEROPTV(-0.8360311073266357942994, 0, 0),  //
                ZEROPTV(-0.6133714327005903973087, 0, 0),  //
                ZEROPTV(-0.3242534234038089290385, 0, 0),  //
                ZEROPTV(0.0, 0, 0),  //
                ZEROPTV(+0.3242534234038089290385, 0, 0),  //
                ZEROPTV(+0.6133714327005903973087, 0, 0),  //
                ZEROPTV(+0.8360311073266357942994, 0, 0),  //
                ZEROPTV(+0.9681602395076260898356, 0, 0)   //
        };
        ///! integration point weights
        static constexpr constexpr_array<double, 9> weights =
                {{
                         0.0812743883615744119719,
                         0.1806481606948574040585,
                         0.2606106964029354623187,
                         0.312347077040002840069,
                         0.330239355001259763165,
                         0.312347077040002840069,
                         0.2606106964029354623187,
                         0.1806481606948574040585,
                         0.0812743883615744119719
                 }};

    };


#endif

}  // namespace TALYFEMLIB