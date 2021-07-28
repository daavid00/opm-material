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
 * \copydoc Opm::DynamicWaParams
 */
#ifndef OPM_DYNAMIC_WA_PARAMS_HPP
#define OPM_DYNAMIC_WA_PARAMS_HPP

#include <opm/material/common/Valgrind.hpp>
#include <opm/material/common/EnsureFinalized.hpp>

#include <cassert>

namespace Opm {

/*!
 * \ingroup FluidMatrixInteractions
 *
 * \brief Specification of the material parameters for the
 *        Dynamic-Wa constitutive relations.
 *
 *\see DynamicWa
 */
template <class TraitsT>
class DynamicWaParams : public EnsureFinalized
{
    typedef typename TraitsT::Scalar Scalar;
public:
    using EnsureFinalized :: finalize;

    typedef TraitsT Traits;

    DynamicWaParams()
    {
        Valgrind::SetUndefined(*this);
    }

    DynamicWaParams(Scalar ePressure, Scalar shapeParam)
        : entryPressure_(ePressure), lambda_(shapeParam)
    {
        finalize();
    }

    /*!
     * \brief Returns the entry pressure [Pa]
     */
    Scalar entryPressure() const
    { EnsureFinalized::check(); return entryPressure_; }

    /*!
     * \brief Set the entry pressure [Pa]
     */
    void setEntryPressure(Scalar v)
    { entryPressure_ = v; }

    /*!
     * \brief Returns the final entry pressure [Pa]
     */
    Scalar finalEntryPressure() const
    { EnsureFinalized::check(); return finalEntryPressure_; }

    /*!
     * \brief Set the final entry pressure [Pa]
     */
    void setFinalEntryPressure(Scalar v)
    { finalEntryPressure_ = v; }

    /*!
     * \brief Returns the beta wa parameter (capillary pressure)
     */
    Scalar beta() const
    { return beta_; }

    /*!
     * \brief Set the beta wa parameter (capillary pressure)
     */
    void setBeta(Scalar v)
    { beta_ = v; }

    /*!
     * \brief Returns the eta wa parameter (relative permeability)
     */
    Scalar eta() const
    { return eta_; }

     /*!
     * \brief Set the eta wa parameter (relative permeability)
     */
    void setEta(Scalar v)
    { eta_ = v; }

    /*!
     * \brief Returns the Ei shape parameter (relative permeability)
     */
    Scalar ei() const
    { return ei_; }

    /*!
    * \brief Set the Ei shape parameter (relative permeability)
    */
    void setEi(Scalar v)
    { ei_ = v; }

    /*!
     * \brief Returns the Ef  shape parameter (relative permeability)
     */
    Scalar ef() const
    { return ef_; }

    /*!
    * \brief Set the Ef shape parameter (relative permeability)
    */
    void setEf(Scalar v)
    { ef_ = v; }

    /*!
     * \brief Returns the lambda shape parameter (capillary pressure)
     */
    Scalar lambda() const
    { EnsureFinalized::check(); return lambda_; }

    /*!
     * \brief Set the lambda shape parameter (capillary pressure)
     */
    void setLambda(Scalar v)
    { lambda_ = v; }

    /*!
     * \brief Returns the Lambda shape parameter (relative permeability)
     */
    Scalar llambda() const
    { EnsureFinalized::check(); return llambda_; }

    /*!
     * \brief Set the Lambda shape parameter (relative permeability)
     */
    void setLlambda(Scalar v)
    { llambda_ = v; }

private:
    Scalar entryPressure_;
    Scalar finalEntryPressure_;
    Scalar lambda_;
    Scalar llambda_;
    Scalar beta_;
    Scalar eta_;
    Scalar ei_;
    Scalar ef_;
};
} // namespace Opm

#endif
