/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFSED.hpp"
#include "FSPSvarIMFSEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* FSPSvarIMFSED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new FSPSvarIMFSEDFamily(this);
}

////////////////////////////////////////////////////////////////////
