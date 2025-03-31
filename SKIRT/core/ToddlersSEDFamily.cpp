/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "ToddlersSEDFamily.hpp"
#include "Constants.hpp"
#include "FatalError.hpp"

////////////////////////////////////////////////////////////////////

ToddlersSEDFamily::ToddlersSEDFamily(SimulationItem* parent, SedMode sedMode, bool includeDust,
                                    Resolution resolution, SFRPeriod sfrPeriod)
{
    parent->addChild(this);
    _sedMode = sedMode;
    _includeDust = includeDust;
    _resolution = resolution;
    _sfrPeriod = sfrPeriod;
    setup();
}

////////////////////////////////////////////////////////////////////

string ToddlersSEDFamily::getResourceNameSuffix() const
{
    string suffix;

    // Add SED mode
    if (_sedMode == SedMode::Cloud)
    {
        suffix += "Cloud_";
    }
    else  // SFRNormalized
    {
        suffix += "SFRNormalized_";
    }

    // Add stellar population parameters based on selected template
    switch (_stellarTemplate)
    {
        case StellarTemplate::SB99Kroupa100Sin:
            suffix += "SB99_kroupa100_sin";
            break;
        case StellarTemplate::BPASSChab100Bin:
            suffix += "BPASS_chab100_bin";
            break;
        case StellarTemplate::BPASSChab300Bin:
            suffix += "BPASS_chab300_bin";
            break;
    }
    suffix += "_";

    // Add dust option
    if (_includeDust)
    {
        suffix += "Dust_";
    }
    else
    {
        suffix += "noDust_";
    }

    // Add resolution
    suffix += _resolution == Resolution::Low ? "lr" : "hr";

    // Add SFR integration period if in SFRNormalized mode
    if (_sedMode == SedMode::SFRNormalized)
    {
        if (_sfrPeriod == SFRPeriod::Period30Myr)
        {
            suffix += "_30Myr";
        }
        else // Period10Myr
        {
            suffix += "_10Myr";
        }
    }

    return suffix;
}

////////////////////////////////////////////////////////////////////

void ToddlersSEDFamily::setupSelfBefore()
{
    SEDFamily::setupSelfBefore();
    
    string name = "ToddlersSEDFamily_";
    name += getResourceNameSuffix();

    if (_sedMode == SedMode::Cloud)
    {
        _cloudTable.open(this, name, "lambda(m),time(Myr),Z(1),SFE(1),n_cl(1/cm3),M_cl(Msun)", "Llambda(W/m)", false);
    }
    else  // SFRNormalized
    {
        _sfrNormalizedTable.open(this, name, "lambda(m),Z(1),SFE(1),n_cl(1/cm3)", "Llambda(W/m)", false);
    }
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> ToddlersSEDFamily::parameterInfo() const
{
    if (_sedMode == SedMode::Cloud)
    {
        return {
            SnapshotParameter::age(),
            SnapshotParameter::metallicity(),
            SnapshotParameter::custom("Star formation efficiency"),
            SnapshotParameter::custom("Cloud number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("Mass", "mass", "Msun"),
            SnapshotParameter::custom("scaling"),
        };
    }
    else  // SFRNormalized
    {
        return {
            SnapshotParameter::metallicity(),
            SnapshotParameter::custom("Star formation efficiency"),
            SnapshotParameter::custom("Cloud number density", "numbervolumedensity", "1/cm3"),
            SnapshotParameter::custom("star formation rate", "massrate", "Msun/yr"),
        };
    }
}

////////////////////////////////////////////////////////////////////

Range ToddlersSEDFamily::intrinsicWavelengthRange() const
{
    if (_sedMode == SedMode::Cloud)
    {
        return _cloudTable.axisRange<0>();
    }
    else  // SFRNormalized
    {
        return _sfrNormalizedTable.axisRange<0>();
    }
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    if (_sedMode == SedMode::Cloud)
    {
        // Extract parameters for Cloud mode
        double age = parameters[0] / (1e6 * Constants::year());  // Convert from s to Myr
        double Z = parameters[1];  // Metallicity
        double SFE = parameters[2];  // Star formation efficiency
        double n_cl = parameters[3] / 1e6;  // Convert from 1/m³ to 1/cm³
        double M_cl = parameters[4] / Constants::Msun();  // Convert from kg to Msun
        double scaling = parameters[5];  // Optional scaling factor

        // Use 6D table: lambda, time, Z, SFE, n_cl, M_cl
        return scaling * _cloudTable(wavelength, age, Z, SFE, n_cl, M_cl);
    }
    else  // SFRNormalized
    {
        // Extract parameters for SFRNormalized mode
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
        double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from Msun/yr to kg/s

        // Use 4D table: lambda, Z, SFE, n_cl
        return sfr * _sfrNormalizedTable(wavelength, Z, SFE, n_cl);
    }
}

////////////////////////////////////////////////////////////////////

double ToddlersSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                             const Array& parameters) const
{
    if (_sedMode == SedMode::Cloud)
    {
        // Extract parameters for Cloud mode
        double age = parameters[0] / (1e6 * Constants::year());  // Convert from s to Myr
        double Z = parameters[1];  // Metallicity
        double SFE = parameters[2];  // Star formation efficiency
        double n_cl = parameters[3] / 1e6;  // Convert from 1/m³ to 1/cm³
        double M_cl = parameters[4] / Constants::Msun();  // Convert from kg to Msun
        double scaling = parameters[5];  // Optional scaling factor

        // Use 6D table: lambda, time, Z, SFE, n_cl, M_cl
        return scaling * _cloudTable.cdf(lambdav, pv, Pv, wavelengthRange, age, Z, SFE, n_cl, M_cl);
    }
    else  // SFRNormalized
    {
        // Extract parameters for SFRNormalized mode
        double Z = parameters[0];
        double SFE = parameters[1];
        double n_cl = parameters[2] / 1e6;  // Convert from 1/m³ to 1/cm³
        double sfr = parameters[3] / Constants::Msun() * Constants::year();  // Convert from Msun/yr to kg/s

        // Use 4D table: lambda, Z, SFE, n_cl
        return sfr * _sfrNormalizedTable.cdf(lambdav, pv, Pv, wavelengthRange, Z, SFE, n_cl);
    }
}

////////////////////////////////////////////////////////////////////