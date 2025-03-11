/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSFRNORMALIZEDSEDFAMILY_HPP
#define TODDLERSSFRNORMALIZEDSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** The ToddlersSFRNormalizedSEDFamily class represents SEDs for star-forming regions,
    normalized by star formation rate. The SEDs are derived from TODDLERS model grid
    calculations. */
class ToddlersSFRNormalizedSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the template and its characteristics */
    ENUM_DEF(TemplateType, SB99Kroupa100Sin, BPASSChab100Bin, BPASSChab300Bin)
        ENUM_VAL(TemplateType, SB99Kroupa100Sin, "Starburst99 with Kroupa IMF (0.1-100 Msun) and single star evolution")
        ENUM_VAL(TemplateType, BPASSChab100Bin, "BPASS with Chabrier IMF (0.1-100 Msun) and binary star evolution")
        ENUM_VAL(TemplateType, BPASSChab300Bin, "BPASS with Chabrier IMF (0.1-300 Msun) and binary star evolution")
    ENUM_END()

    /** The enumeration type indicating the presence of dust */
    ENUM_DEF(Dust, Yes, No)
        ENUM_VAL(Dust, Yes, "Dust is present in SF regions")
        ENUM_VAL(Dust, No, "No dust is present in SF regions, uses incident stellar continuum")
    ENUM_END()

    /** The enumeration type indicating the wavelength resolution */
    ENUM_DEF(Resolution, Low, High)
        ENUM_VAL(Resolution, Low, "Low wavelength resolution (continuum and lines at R=300)")
        ENUM_VAL(Resolution, High, "High wavelength resolution (continuum at R=300 and lines at R=5e4)")
    ENUM_END()

    /** The enumeration type indicating the SFR integration period */
    ENUM_DEF(SFRPeriod, Period10Myr, Period30Myr)
        ENUM_VAL(SFRPeriod, Period10Myr, "SFR integrated over 10 Myr (default)")
        ENUM_VAL(SFRPeriod, Period30Myr, "SFR integrated over 30 Myr")
    ENUM_END()

    ITEM_CONCRETE(ToddlersSFRNormalizedSEDFamily, SEDFamily,
                 "a TODDLERS SFR-normalized SED family for emission from star-forming regions")
        PROPERTY_ENUM(templateType, TemplateType, "the stellar template, IMF, and evolution model to use")
        ATTRIBUTE_DEFAULT_VALUE(templateType, "SB99Kroupa100Sin")

        PROPERTY_ENUM(dust, Dust, "the presence of dust")
        ATTRIBUTE_DEFAULT_VALUE(dust, "Yes")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

        PROPERTY_ENUM(sfrPeriod, SFRPeriod, "the SFR integration time period")
        ATTRIBUTE_DEFAULT_VALUE(sfrPeriod, "Period10Myr")
    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family. The newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy, and its setup() function has been called. */
    explicit ToddlersSFRNormalizedSEDFamily(SimulationItem* parent, Dust dust, Resolution resolution, 
                                          SFRPeriod sfrPeriod = SFRPeriod::Period10Myr);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each object specifies unit information
        and a human-readable description for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family from the stored table. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (radiative power per
        unit of wavelength) for the %SED with the specified parameters at the specified wavelength. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs the normalized probability density function (pdf) and
        cumulative distribution function (cdf) for the %SED with the specified parameters. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
              const Array& parameters) const override;

private:
    /** Returns the filename suffix for the current configuration. */
    string getResourceNameSuffix() const;

    //======================== Data Members =======================

private:
    StoredTable<4> _table;  // 4D table: wavelength, Z, SFE, n_cl
};

//////////////////////////////////////////////////////////////////////

#endif
