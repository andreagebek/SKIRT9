/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef BPASSCHABRIER100SED_HPP
#define BPASSCHABRIER100SED_HPP

#include "FamilySED.hpp"

////////////////////////////////////////////////////////////////////

/** BpassChabrier100SED is a class that represents spectral energy distributions of simple stellar populations
    (SSPs), parameterized on metallicity and age according to the BPASS model that includes binary
    stellar systems and that assumes a Chabrier IMF with a lower mass limit of 1 solar mass and an upper mass
    limit of 100 solar masses. See the BpassChabrier100SEDFamily class for more information. */

class BpassChabrier100SED : public FamilySED
{
    ITEM_CONCRETE(BpassChabrier100SED, FamilySED, "a BPASS single stellar population SED with an upper mass limit of 100 solar masses")

        PROPERTY_DOUBLE(metallicity, "the metallicity of the SSP")
        ATTRIBUTE_MIN_VALUE(metallicity, "[1e-5")
        ATTRIBUTE_MAX_VALUE(metallicity, "0.04]")
        ATTRIBUTE_DEFAULT_VALUE(metallicity, "0.02")

        PROPERTY_DOUBLE(age, "the age of the SSP")
        ATTRIBUTE_QUANTITY(age, "time")
        ATTRIBUTE_MIN_VALUE(age, "[1 Myr")
        ATTRIBUTE_MAX_VALUE(age, "100 Gyr]")
        ATTRIBUTE_DEFAULT_VALUE(age, "5 Gyr")

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
