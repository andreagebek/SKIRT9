/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProGenySED.hpp"
#include "ProGenySEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* ProGenySED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age, _alpha);

    // construct the library of SED models
    return new ProGenySEDFamily(this);
}

////////////////////////////////////////////////////////////////////
