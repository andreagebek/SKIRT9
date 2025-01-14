/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

FSPSvarIMFSEDFamily::FSPSvarIMFSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void FSPSvarIMFSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "FSPSSEDFamily_Variable", "lambda(m),Z(1),alpha(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FSPSvarIMFSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::custom("IMF slope"), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range FSPSvarIMFSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];
    double t = parameters[3] / Constants::year();

    return M * _table(wavelength, Z, alpha, t);
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                          const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];
    double t = parameters[3] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, alpha, t);
}

////////////////////////////////////////////////////////////////////
