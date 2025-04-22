/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "BpassChabrier100SED.hpp"
#include "BpassChabrier100SEDFamily.hpp"
#include "NR.hpp"

////////////////////////////////////////////////////////////////////

const SEDFamily* BpassChabrier100SED::getFamilyAndParameters(Array& parameters)
{
    // set the parameters using arbitrary scaling
    NR::assign(parameters, 1., _metallicity, _age);

    // construct the library of SED models
    return new BpassChabrier100SEDFamily(this);
}

////////////////////////////////////////////////////////////////////
