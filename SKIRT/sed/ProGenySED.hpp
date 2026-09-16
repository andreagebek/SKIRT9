/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROGENYSED_HPP
#define PROGENYSED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** A ProGenySED class instance represents a simple stellar population %SED taken from the
    ProGeny code (Robotham & Bellstedt 2024) with a Chabrier IMF (0.1-100 \f$\mathrm{M}_\odot\f$).
    This particular variant is based on the sMILES spectral library which features different alpha-
    enhancements (PG_Ch_Mi_Kn_Alpha_X), shared by Themiya Nanayakkara in 2026.
    The SED is parametrized on stellar metallicity, age, and alpha-enhancement.
    See the ProGenySEDFamily class for more information. */
class ProGenySED : public FamilySED
{
    ITEM_CONCRETE(ProGenySED, FamilySED, "a ProGeny simple stellar population SED")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[1e-6")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.06]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[0.1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "20 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

        PROPERTY_DOUBLE(alpha, "the alpha-enhancement of the SSP")
        ATTRIBUTE_MIN_VALUE(alpha, "-0.2")
        ATTRIBUTE_MAX_VALUE(alpha, "0.6")
        ATTRIBUTE_DEFAULT_VALUE(alpha, "0")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function returns a newly created SEDFamily object (which is already hooked into the
        simulation item hierachy so it will be automatically deleted) and stores the parameters for
        the specific %SED configured by the user in the specified array. */
    const SEDFamily* getFamilyAndParameters(Array& parameters) override;
};

////////////////////////////////////////////////////////////////////

#endif
