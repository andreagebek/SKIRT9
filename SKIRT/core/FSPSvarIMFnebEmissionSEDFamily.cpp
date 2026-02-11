/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSvarIMFnebEmissionSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

FSPSvarIMFnebEmissionSEDFamily::FSPSvarIMFnebEmissionSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void FSPSvarIMFnebEmissionSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "FSPSSEDFamily_Variable_withNebularEmission", "lambda(m),Z(1),alpha(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FSPSvarIMFnebEmissionSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::custom("IMF slope"), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range FSPSvarIMFnebEmissionSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFnebEmissionSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];
    double t = parameters[3] / Constants::year();

    return M * _table(wavelength, Z, alpha, t);
}

////////////////////////////////////////////////////////////////////

double FSPSvarIMFnebEmissionSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                          const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double alpha = parameters[2];
    double t = parameters[3] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, alpha, t);
}

////////////////////////////////////////////////////////////////////
