/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "IsotropicAngularDistribution.hpp"
#include "Random.hpp"

////////////////////////////////////////////////////////////////////

int IsotropicAngularDistribution::dimension() const
{
    return 1;
}

////////////////////////////////////////////////////////////////////

double IsotropicAngularDistribution::probabilityForDirection(Direction /*bfk*/) const
{
    return 1.;
}

//////////////////////////////////////////////////////////////////////

Direction IsotropicAngularDistribution::generateDirection() const
{
    return random()->direction();
}

////////////////////////////////////////////////////////////////////
