/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSFRNormalizedSEDFamily.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

ToddlersSFRNormalizedSEDFamily::ToddlersSFRNormalizedSEDFamily(SimulationItem* parent, Dust dust,
                                                              Resolution resolution) {
    parent->addChild(this);
    _dust = dust;
    _resolution = resolution;
    setup();
}

////////////////////////////////////////////////////////////////////

void ToddlersSFRNormalizedSEDFamily::validateConfiguration() const 
{
    // SB99 is only valid with kroupa100 and single stars
    if (_stellarTemplate == StellarTemplate::SB99)
    {
        if (_imf != IMF::kroupa100 || _starType != StarType::sin)
        { 
            throw  FATALERROR("SB99 models are only available with Kroupa IMF and single star evolution");
        }
    }
    
    // BPASS is only valid with Chabrier IMFs and binary stars
    if (_stellarTemplate == StellarTemplate::BPASS)
    {
        if (_starType != StarType::bin)
        {
            throw FATALERROR("BPASS models are only available with binary star evolution");
        }
        if (_imf != IMF::chab100 && _imf != IMF::chab300)
        {
            throw FATALERROR("BPASS models are only available with Chabrier IMF");
        }
    }
}

////////////////////////////////////////////////////////////////////

string ToddlersSFRNormalizedSEDFamily::getResourceNameSuffix() const
{
    validateConfiguration();
    
    string suffix;

    // Add stellar population parameters
    switch (_stellarTemplate)
    {
        case StellarTemplate::SB99: suffix += "SB99"; break;
        case StellarTemplate::BPASS: suffix += "BPASS"; break;
    }
    suffix += "_";

    switch (_imf)
    {
        case IMF::kroupa100: suffix += "kroupa100"; break;
        case IMF::chab100: suffix += "chab100"; break;
        case IMF::chab300: suffix += "chab300"; break;
    }
    suffix += "_";

    switch (_starType)
    {
        case StarType::sin: suffix += "sin"; break;
        case StarType::bin: suffix += "bin"; break;
    }
    suffix += "_";

    // Add dust option
    if (_dust == Dust::No)
    {
        suffix += "noDust_";
    }

    // Add resolution
    suffix += _resolution == Resolution::Low ? "lr" : "hr";

    return suffix;
}

////////////////////////////////////////////////////////////////////

void ToddlersSFRNormalizedSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();
    
    validateConfiguration();
    
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
    double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from kg/s to Msun/yr

    return sfr * _table.cdf(lambdav, pv, Pv, wavelengthRange, Z, SFE, n_cl);
}

///////////////////////////////////////////////////////////////////