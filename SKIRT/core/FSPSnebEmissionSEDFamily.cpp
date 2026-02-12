/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "FSPSnebEmissionSEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

FSPSnebEmissionSEDFamily::FSPSnebEmissionSEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void FSPSnebEmissionSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "FSPSSEDFamily_Chabrier_withNebularEmission", "lambda(m),Z(1),t(yr)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> FSPSnebEmissionSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age()};
}

////////////////////////////////////////////////////////////////////

Range FSPSnebEmissionSEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double FSPSnebEmissionSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table(wavelength, Z, t);
}

////////////////////////////////////////////////////////////////////

double FSPSnebEmissionSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                          const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t);
}

////////////////////////////////////////////////////////////////////
