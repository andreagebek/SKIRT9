/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFcontSFHSED.hpp"
#include "FSPSvarIMFcontSFHSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* FSPSvarIMFcontSFHSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _alpha);

    // construct the library of SED models
    return new FSPSvarIMFcontSFHSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
