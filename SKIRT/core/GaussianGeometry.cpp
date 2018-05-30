/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "GaussianGeometry.hpp"
#include "Random.hpp"

//////////////////////////////////////////////////////////////////////

void GaussianGeometry::setupSelfBefore()
{
    SepAxGeometry::setupSelfBefore();

    // calculate cached values
    _rho0 = 1.0 / pow(sqrt(2.0*M_PI)*_sigma,3) / _q;
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::density(double R, double z) const
{
    double m2 = R*R + z*z/(_q*_q);
    double sigma2 = _sigma*_sigma;
    return _rho0 * exp(-0.5*m2/sigma2);
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::randomCylRadius() const
{
    double X = random()->uniform();
    return _sigma * sqrt(-2.0*log(X));
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::randomZ() const
{
    return (_q*_sigma) * random()->gauss();
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::SigmaR() const
{
    return 1.0/(4.0*M_PI*_q*_sigma*_sigma);
}

////////////////////////////////////////////////////////////////////

double GaussianGeometry::SigmaZ() const
{
    return 1.0/(2.0*M_PI*_sigma*_sigma);
}

////////////////////////////////////////////////////////////////////