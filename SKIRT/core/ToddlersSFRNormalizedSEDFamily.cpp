/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSFRNormalizedSEDFamily.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

ToddlersSFRNormalizedSEDFamily::ToddlersSFRNormalizedSEDFamily(SimulationItem* parent, Dust dust,
                                                              Resolution resolution, SFRPeriod sfrPeriod) {
    parent->addChild(this);
    _dust = dust;
    _resolution = resolution;
    _sfrPeriod = sfrPeriod;
    setup();
}

////////////////////////////////////////////////////////////////////

string ToddlersSFRNormalizedSEDFamily::getResourceNameSuffix() const
{
    string suffix;

    // Add stellar population parameters based on selected template
    switch (_templateType)
    {
        case TemplateType::SB99Kroupa100Sin:
            suffix += "SB99_kroupa100_sin";
            break;
        case TemplateType::BPASSChab100Bin:
            suffix += "BPASS_chab100_bin";
            break;
        case TemplateType::BPASSChab300Bin:
            suffix += "BPASS_chab300_bin";
            break;
    }
    suffix += "_";

    // Add dust option
    if (_dust == Dust::No)
    {
        suffix += "noDust_";
    }

    // Add resolution
    suffix += _resolution == Resolution::Low ? "lr" : "hr";

    // Add SFR integration period
    suffix += _sfrPeriod == SFRPeriod::Period10Myr ? "_10Myr" : "_30Myr";

    return suffix;
}

////////////////////////////////////////////////////////////////////

void ToddlersSFRNormalizedSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();
    
    string name = "ToddlersSFRNormalizedSEDFamily_";
    name += getResourceNameSuffix();

    _table.open(this, name, "lambda(m),Z(1),SFE(1),n_cl(1/cm3)", "Llambda(W/m)", false);
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ToddlersSFRNormalizedSEDFamily::parameterInfo() const {
    return {
        SnapshotParameter::metallicity(),
        SnapshotParameter::custom("Star formation efficiency"),
        SnapshotParameter::custom("Cloud number density", "numbervolumedensity", "1/cm3"),
        SnapshotParameter::custom("star formation rate", "massrate", "Msun/yr"),
    };
}

////////////////////////////////////////////////////////////////////

Range ToddlersSFRNormalizedSEDFamily::intrinsicWavelengthRange() const {
    return _table.axisRange<0>();
}

////////////////////////////////////////////////////////////////////

double ToddlersSFRNormalizedSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const {
    double Z = parameters[0];
    double SFE = parameters[1];
    double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
    double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from Msun/yr to kg/s

    return sfr * _table(wavelength, Z, SFE, n_cl);
}

////////////////////////////////////////////////////////////////////

double ToddlersSFRNormalizedSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                                          const Array& parameters) const {
    double Z = parameters[0];
    double SFE = parameters[1];
    double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
    double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from Msun/yr to kg/s

    return sfr * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, SFE, n_cl);
}

///////////////////////////////////////////////////////////////////
