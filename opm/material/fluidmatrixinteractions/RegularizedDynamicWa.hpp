// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 * \copydoc Opm::RegularizedDynamicWa
 */
#ifndef REGULARIZED_DYNAMIC_WA_HPP
#define REGULARIZED_DYNAMIC_WA_HPP

#include "DynamicWa.hpp"
#include "RegularizedDynamicWaParams.hpp"

#include <opm/material/common/Spline.hpp>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 * \brief Implementation of the regularized Dynamic-Wa capillary
 *        pressure / relative permeability <-> saturation relation.
 *
 * This class bundles the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 * In order to avoid very steep gradients the marginal values are
 * "regularized".  This means that in stead of following the curve of
 * the material law in these regions, some linear approximation is
 * used.  Doing this is not worse than following the material
 * law. E.g. for very low wetting phase values the material laws
 * predict infinite values for \f$p_c\f$ which is completely
 * unphysical. In case of very high wetting phase saturations the
 * difference between regularized and "pure" material law is not big.
 *
 * Regularizing has the additional benefit of being numerically
 * friendly: Newton's method does not like infinite gradients.
 *
 * The implementation is accomplished as follows:
 * - check whether we are in the range of regularization
 *   - yes: use the regularization
 *   - no: forward to the standard material law.
 *
 * \see DynamicWa
 */
template <class TraitsT, class ParamsT = RegularizedDynamicWaParams<TraitsT> >
class RegularizedDynamicWa : public TraitsT
{
    typedef Opm::DynamicWa<TraitsT, ParamsT> DynamicWa;

public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The regularized Dynamic-Wa capillary pressure law only "
                  "applies to the case of two fluid phases");

    //! Specify whether this material law implements the two-phase
    //! convenience API
    static const bool implementsTwoPhaseApi = true;

    //! Specify whether this material law implements the two-phase
    //! convenience API which only depends on the phase saturations
    static const bool implementsTwoPhaseSatApi = true;

    //! Specify whether the quantities defined by this material law
    //! are saturation dependent
    static const bool isSaturationDependent = true;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the absolute pressure
    static const bool isPressureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are temperature dependent
    static const bool isTemperatureDependent = false;

    //! Specify whether the quantities defined by this material law
    //! are dependent on the phase composition
    static const bool isCompositionDependent = false;

    static_assert(Traits::numPhases == 2,
                  "The number of fluid phases must be two if you want to use "
                  "this material law!");

    /*!
     * \brief The capillary pressure-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative pressure of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the capillary pressure
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The relative permeability-saturation curves depending on absolute saturations.
     *
     * \param values A random access container which stores the
     *               relative permeability of each fluid phase.
     * \param params The parameter object expressing the coefficients
     *               required by the material law.
     * \param fs The fluid state for which the relative permeabilities
     *           ought to be calculated
     */
    template <class Container, class FluidState>
    static void relativePermeabilities(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = krw<FluidState, Evaluation>(params, fs);
        values[Traits::nonWettingPhaseIdx] = krn<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief A regularized Dynamic-Wa capillary pressure-saturation curve.
     *
     * This is a regularized variant of the Dynamic-Wa curve:
     *
     * - For wetting phase saturations lower than the threshold saturation, the
     *   \f$p_c(S_w)\f$ curve is extrapolated using a straight line exhibiting the slope
     *   unregularized capillary pressure curve at the threshold saturation.
     * - For wetting phase saturations larger than 1, the curve is extrapolated using a
     *   straight line that exhibits the slope of the unregularized curve at \f$S_w =
     *   1\f$
     *
     * \sa DynamicWa::pcnw
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const auto& Sw = Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        const auto& Wa = Opm::decay<Evaluation>(fs.wa());
        return twoPhaseSatPcnw(params, Sw, Wa);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw, const Evaluation& Wa)
    {
        const Scalar Sthres = params.pcnwLowSw();

        if (Sw <= Sthres) {
            Scalar pcnw_SwLow = DynamicWa::twoPhaseSatPcnw(params, Sthres, Wa);
            const Scalar eps = 1e-7;
            Scalar delta = 0.0;
            Scalar pc1;
            Scalar pc2;
            pc2 = DynamicWa::twoPhaseSatPcnw(params, Sthres + eps, Wa);
            delta += eps;
            pc1 = DynamicWa::twoPhaseSatPcnw(params, Sthres - eps, Wa);
            delta += eps;
            Scalar m = (pc2 - pc1)/delta;

            return pcnw_SwLow + m*(Sw - Sthres);
        }
        else if (Sw >= 1.0) {
            Scalar pcnw_SwHigh = DynamicWa::twoPhaseSatPcnw(params, 1.0, Wa);
            const Scalar eps = 1e-7;
            Scalar delta = 0.0;
            Scalar pc1;
            Scalar pc2;
            pc2 = DynamicWa::twoPhaseSatPcnw(params, 1.0, Wa);
            pc1 = DynamicWa::twoPhaseSatPcnw(params, 1.0 - eps, Wa);
            delta += eps;
            Scalar m = (pc2 - pc1)/delta;
            return pcnw_SwHigh + m*(Sw - 1.0);
        }

        // if the effective saturation is in an 'reasonable'
        // range, we use the real Dynamic-Wa saturation functions...
        return DynamicWa::twoPhaseSatPcnw(params, Sw, Wa);
    }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        wetting phase of the Dynamic-Wa curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Dynamic-Wa wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 0
     * - For wetting phase saturations larger than 1, return 1
     *
     * \sa DynamicWa::krw
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const auto& Sw = Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        const auto& Wa = Opm::decay<Evaluation>(fs.wa());
        return twoPhaseSatKrw(params, Sw, Wa);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw, const Evaluation& Wa)
    {
        if (Sw <= 0.0)
            return 0.0;
        else if (Sw >= 1.0)
            return 1.0;

        return DynamicWa::twoPhaseSatKrw(params, Sw, Wa);
    }

    /*!
     * \brief Regularized version of the relative permeability of the
     *        non-wetting phase of the Dynamic-Wa curves.
     *
     * The approach for regularization is very similar to the one of
     * the capillary pressure, but it does not avoid kinks:
     * - For wetting phase saturations between 0 and 1, use the
     *   unregularized Dynamic-Wa non-wetting phase relative
     *   permeability
     * - For wetting phase saturations smaller than 0, return 1
     * - For wetting phase saturations larger than 1, return 0
     *
     * \sa DynamicWa::krn
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krn(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            1.0 - Opm::decay<Evaluation>(fs.saturation(Traits::nonWettingPhaseIdx));
        const auto& Wa = Opm::decay<Evaluation>(fs.wa());
        return twoPhaseSatKrn(params, Sw, Wa);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrn(const Params& params, const Evaluation& Sw, const Evaluation& Wa)
    {
        if (Sw >= 1.0)
            return 0.0;
        else if (Sw <= 0.0)
            return 1.0;

        return DynamicWa::twoPhaseSatKrn(params, Sw, Wa);
    }
};
} // namespace Opm

#endif
