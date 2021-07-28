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
 * \copydoc Opm::DynamicWa
 */
#ifndef OPM_DYNAMIC_WA_HPP
#define OPM_DYNAMIC_WA_HPP

#include "DynamicWaParams.hpp"

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Exceptions.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>

namespace Opm {
/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Implementation of the Dynamic-Wa capillary pressure <->
 *        saturation relation.
 *
 * This class provides the "raw" curves as static members and doesn't
 * concern itself converting absolute to effective saturations and
 * vice versa.
 *
 *\see DynamicWaParams
 */
template <class TraitsT, class ParamsT = DynamicWaParams<TraitsT> >
class DynamicWa : public TraitsT
{
public:
    typedef TraitsT Traits;
    typedef ParamsT Params;
    typedef typename Traits::Scalar Scalar;

    //! The number of fluid phases to which this material law applies.
    static const int numPhases = Traits::numPhases;
    static_assert(numPhases == 2,
                  "The Dynamic-Wa capillary pressure law only applies "
                  "to the case of two fluid phases");

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
     * \brief The capillary pressure-saturation curves.
     */
    template <class Container, class FluidState>
    static void capillaryPressures(Container& values, const Params& params, const FluidState& fs)
    {
        typedef typename std::remove_reference<decltype(values[0])>::type Evaluation;

        values[Traits::wettingPhaseIdx] = 0.0; // reference phase
        values[Traits::nonWettingPhaseIdx] = pcnw<FluidState, Evaluation>(params, fs);
    }

    /*!
     * \brief The relative permeability-saturation curves.
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
     * \brief The capillary pressure-saturation curve in the
     *        Dynamic-Wa model ().
     *
     * The empirical Dynamic-Wa capillary pressure-saturation
     * function is defined as
     * \f[
     * p_C = p_e\overline{S}_w^{-1/\lambda}
     * \f]
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Dynamic-Wa: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation pcnw(const Params& params, const FluidState& fs)
    {
        const Evaluation& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        const Evaluation& Wa =
            Opm::decay<Evaluation>(fs.wa());

        assert(0.0 <= Sw && Sw <= 1.0);

        return twoPhaseSatPcnw(params, Sw, Wa);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatPcnw(const Params& params, const Evaluation& Sw, const Evaluation& Wa)
    {
        assert(0.0 <= Sw && Sw <= 1.0);

        return (1.0 + (params.finalEntryPressure()/params.entryPressure()-1.0) * (Sw*Wa) / (params.beta() + Sw*Wa)) * params.entryPressure()*Opm::pow(Sw, -1/params.lambda());
    }

    /*!
     * \brief The relative permeability for the wetting phase of
     *        the medium implied by the Dynamic-Wa
     *        parameterization.
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Dynamic-Wa: Entry pressure and shape factor)
     */
    template <class FluidState, class Evaluation = typename FluidState::Scalar>
    static Evaluation krw(const Params& params, const FluidState& fs)
    {
        const auto& Sw =
            Opm::decay<Evaluation>(fs.saturation(Traits::wettingPhaseIdx));
        const auto& Wa = Opm::decay<Evaluation>(fs.wa());

        return twoPhaseSatKrw(params, Sw, Wa);
    }

    template <class Evaluation>
    static Evaluation twoPhaseSatKrw(const Params& params, const Evaluation& Sw, const Evaluation& Wa)
    {
        assert(0.0 <= Sw && Sw <= 1.0);

        return Opm::min(params.eta()*(Wa) + params.ei(),params.ef()) * Opm::pow(Sw, params.llambda())/(1.0- Sw + Opm::min(params.eta()*(Wa) + params.ei(),params.ef()) * Opm::pow(Sw, params.llambda()));
    }

    /*!
     * \brief The relative permeability for the non-wetting phase of
     *        the medium as implied by the Dynamic-Wa
     *        parameterization.
     *
     * \param params The parameters of the capillary pressure curve
     *               (for Dynamic-Wa: Entry pressure and shape factor)
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
        assert(0.0 <= Sw && Sw <= 1.0);
        const Evaluation Sn = 1.0 - Sw;

        return (1 - Sw) / (1.0 - Sw + Opm::min(params.eta()*(Wa) + params.ei() , params.ef()) * Opm::pow(Sw, params.llambda()));
    }

};
} // namespace Opm

#endif // DYNAMIC_WA_HPP
