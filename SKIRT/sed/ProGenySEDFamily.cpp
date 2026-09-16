/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ProGenySEDFamily.hpp"
#include "Constants.hpp"

////////////////////////////////////////////////////////////////////

ProGenySEDFamily::ProGenySEDFamily(SimulationItem* parent)
{
    parent->addChild(this);
    setup();
}

////////////////////////////////////////////////////////////////////

void ProGenySEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();

    _table.open(this, "ProGenySEDFamily", "lambda(m),Z(1),t(yr),alpha(1)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ProGenySEDFamily::parameterInfo() const
{
    return {SnapshotParameter::initialMass(), SnapshotParameter::metallicity(), SnapshotParameter::age(), SnapshotParameter::custom("alpha-enhancement")};
}

////////////////////////////////////////////////////////////////////

Range ProGenySEDFamily::intrinsicWavelengthRange() const
{
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double ProGenySEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();
    double alpha = parameters[3];

    return M * _table(wavelength, Z, t, alpha);
}

////////////////////////////////////////////////////////////////////

double ProGenySEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                 const Array& parameters) const
{
    double M = parameters[0] / Constants::Msun();
    double Z = parameters[1];
    double t = parameters[2] / Constants::year();
    double alpha = parameters[3];

    return M * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, t, alpha);
}

////////////////////////////////////////////////////////////////////
