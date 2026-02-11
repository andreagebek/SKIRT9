/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFnebEmissionSED.hpp"
#include "FSPSvarIMFnebEmissionSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* FSPSvarIMFnebEmissionSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _alpha, _age);

    // construct the library of SED models
    return new FSPSvarIMFnebEmissionSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
