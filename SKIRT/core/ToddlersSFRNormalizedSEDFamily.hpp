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
    /** The enumeration type indicating the stellar population synthesis model */
    ENUM_DEF(StellarTemplate, SB99, BPASS)
        ENUM_VAL(StellarTemplate, SB99, "Uses Starburst99 stellar population models")
        ENUM_VAL(StellarTemplate, BPASS, "Uses BPASS stellar population models")
    ENUM_END()

    /** The enumeration type indicating the initial mass function */
    ENUM_DEF(IMF, kroupa100, chab100, chab300)
        ENUM_VAL(IMF, kroupa100, "Kroupa IMF from 0.1 to 100 Msun")
        ENUM_VAL(IMF, chab100, "Chabrier IMF from 0.1 to 100 Msun")
        ENUM_VAL(IMF, chab300, "Chabrier IMF from 0.1 to 300 Msun")
    ENUM_END()

    /** The enumeration type indicating the stellar population type */
    ENUM_DEF(StarType, sin, bin)
        ENUM_VAL(StarType, sin, "Single star evolution")
        ENUM_VAL(StarType, bin, "Binary star evolution")
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

    ITEM_CONCRETE(ToddlersSFRNormalizedSEDFamily, SEDFamily,
                 "a TODDLERS SFR-normalized SED family for emission from star-forming regions")
        PROPERTY_ENUM(stellarTemplate, StellarTemplate, "the stellar population synthesis model")
        ATTRIBUTE_DEFAULT_VALUE(stellarTemplate, "SB99")

        PROPERTY_ENUM(imf, IMF, "the initial mass function")
        ATTRIBUTE_DEFAULT_VALUE(imf, "kroupa100")

        PROPERTY_ENUM(starType, StarType, "the stellar population type")
        ATTRIBUTE_DEFAULT_VALUE(starType, "sin")

        PROPERTY_ENUM(dust, Dust, "the presence of dust")
        ATTRIBUTE_DEFAULT_VALUE(dust, "Yes")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")
    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family. The newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy, and its setup() function has been called. */
    explicit ToddlersSFRNormalizedSEDFamily(SimulationItem* parent, Dust dust, Resolution resolution);

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

    /** Validates that the template/IMF/star type combination is valid */
    void validateConfiguration() const;
};

//////////////////////////////////////////////////////////////////////

#endif