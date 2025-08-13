/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFcontSFHSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

FSPSvarIMFcontSFHSEDFamily::FSPSvarIMFcontSFHSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void FSPSvarIMFcontSFHSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "FSPSSEDFamily_Variable_CSFH10Myr", "lambda(m),Z(1),alpha(1)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FSPSvarIMFcontSFHSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("formed mass", "mass", "Msun"), SnapshotParameter::metallicity(), SnapshotParameter::custom("IMF slope")};
}

////////////////////////////////////////////////////////////////////

Range FSPSvarIMFcontSFHSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFcontSFHSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];

    return M * _table(wavelength, Z, alpha);
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFcontSFHSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                          const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, alpha);
}

////////////////////////////////////////////////////////////////////
