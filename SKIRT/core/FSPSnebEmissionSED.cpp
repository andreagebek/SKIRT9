/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSnebEmissionSED.hpp"
#include "FSPSnebEmissionSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* FSPSnebEmissionSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new FSPSnebEmissionSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
