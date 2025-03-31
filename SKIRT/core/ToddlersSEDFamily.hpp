/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef TODDLERSSEDFAMILY_HPP
#define TODDLERSSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ToddlersSEDFamily class represents the family of star-forming 
    region templates from the TODDLERS (Time evolution of Observables including Dust Diagnostics 
    and Line Emission from Regions containing young Stars) model suite.
    
    The TODDLERS model calculates the spherical evolution of a gas cloud around a young stellar 
    cluster, accounting for stellar feedback processes such as stellar winds, supernovae, radiation pressure, 
    and gravitational forces (see https://ui.adsabs.harvard.edu/abs/2023MNRAS.526.3871K and
    https://ui.adsabs.harvard.edu/abs/2024A&A...692A..79K). 

    This class offers several configuration options through its properties:
    
    - sedMode: Determines how the SED is calculated and scaled:
      - Cloud: SEDs for individual star-forming clouds with explicit time evolution
      - SFRNormalized: SEDs pre-integrated over time and cloud mass spectrum, directly scaled by SFR
    
    - stellarTemplate: Determines the stellar population model, IMF, and stellar evolution:
      - SB99Kroupa100Sin: Starburst99 models with Kroupa IMF (0.1-100 \f$\mathrm{M}_\odot\f$) and single star evolution
      - BPASSChab100Bin: BPASS models with Chabrier IMF (0.1-100 \f$\mathrm{M}_\odot\f$) and binary star evolution
      - BPASSChab300Bin: BPASS models with Chabrier IMF (0.1-300 \f$\mathrm{M}_\odot\f$) and binary star evolution
    
    - includeDust: Boolean flag determining dust processing in the SEDs:
      - true (default): 
        - In low resolution: Uses the total emission reported by Cloudy, which includes attenuated stellar, nebular,
         and dust continuum components.
        - In high resolution: Uses dust-attenuated emission lines (emergent luminosities) added to 
          the total emission after removing the low-resolution lines. The high-resolution mode includes a limited set of 
          approximately 140 emission lines tracked by TODDLERS, so the replacement is not one-to-one. Emergent luminosity values are calculated 
          using escape probabilities for diffuse radiation in Cloudy. Roughly, these drop off as \f$E_2(\tau)\f$, where \f$E_2\f$ 
          is the exponential integral and \f$\tau\f$ is the optical depth. \f$E_2(\tau)\f$ drops off faster than \f$\exp(-\tau)\f$. 
          This is discussed in section 3.1.1 of https://ui.adsabs.harvard.edu/abs/2007A%26A...467..187R.
      - false: 
        - In low resolution: Uses only the incident stellar radiation field (stellar continuum) without 
          any gas or dust processing.
        - In high resolution: Uses intrinsic emission line luminosities (without foreground attenuation) added to 
          the incident stellar continuum.
    
    - resolution: Spectral resolution of the SEDs:
      - Low (default): Entire spectrum (continuum and lines) at \f$R=300\f$ resolution
      - High: Low-resolution continuum with selected emission lines represented as high-resolution 
        Gaussian profiles (\f$R = \lambda/\Delta\lambda_G = 5\times 10^4\f$) sampled using 37 points 
        per line. This approach offers better line-to-continuum contrast while maintaining computational efficiency
        when running a large set of Cloudy models.
    
    - sfrPeriod: Time period over which SFR is averaged and integrated (only used in SFRNormalized mode):
    - Period10Myr: SFR integrated over 10 Myr (default)
    - Period30Myr: SFR integrated over 30 Myr
    Note: This property will appear in the ski file even when in Cloud mode, but it has no effect in that mode.
    The value is only relevant and used when sedMode is set to SFRNormalized.

    The SEDs are pre-computed for different combinations of:
    1. Evolution time (Cloud mode only): The time since the start of evolution
       - Time-dependent evolution from 0.1 to 30 Myr
       - Note: At a given evolution time, multiple stellar populations of different ages 
         may be present if recollapse has occurred, triggering subsequent generations of star formation
    2. Metallicity (Z): The parameter space depends on the stellar template:
       - SB99: Range from Z=0.001 to Z=0.04 (5 values)
       - BPASS: Range from Z=0.001 to Z=0.04 (11 values)
    3. Star formation efficiency (SFE): The fraction of cloud mass converted to stars
       - SB99: Range from 0.01 to 0.15 (7 values)
       - BPASS: Range from 0.01 to 0.1 (5 values)
    4. Cloud number density: The initial density of the star-forming cloud
       - SB99: Range from 10 to 2560 \f$\mathrm{cm}^{-3}\f$ (9 values)
       - BPASS: Range from 40 to 640 \f$\mathrm{cm}^{-3}\f$ (5 values)
    5. Cloud mass (Cloud mode only): The mass of the star-forming cloud
       - Range from 10^5 to 10^6.75 \f$\mathrm{M}_\odot\f$

    When using Cloud mode without recollapse effects (e.g., at early times or with parameters 
    that don't trigger recollapse), the resulting SEDs match the shape of standard mass-scaled 
    stellar templates already available in SKIRT. However, when recollapse occurs, the presence 
    of multiple stellar populations of different ages creates a distinct SED shape.

    When using Cloud mode, the parameters must appear in the following order, with the 
    specified default units unless overridden by column header info:

    \f[ t\,(\mathrm{Myr}) \quad Z\,(\mathrm{dimensionless}) \quad 
    \mathrm{SFE}\,(\mathrm{dimensionless}) \quad n_{\text{cl}}\,(\mathrm{cm}^{-3}) \quad 
    M_{\text{cl}}\,(\mathrm{M}_\odot) \quad \mathrm{scaling}\,(\mathrm{dimensionless}) \f]

    where \f$t\f$ is the evolution time since the start of star formation, \f$Z\f$ is the metallicity, \f$\mathrm{SFE}\f$ is the 
    star formation efficiency, \f$n_{\text{cl}}\f$ is the cloud number density, \f$M_{\text{cl}}\f$ is the cloud mass,
    and \f$\mathrm{scaling}\f$ is an arbitrary scaling factor (typically 1).
    
    When using SFRNormalized mode, the following parameters must appear in order:

    \f[ Z\,(\mathrm{dimensionless}) \quad \mathrm{SFE}\,(\mathrm{dimensionless}) \quad 
    n_{\text{cl}}\,(\mathrm{cm}^{-3}) \quad \dot{M}_*\,(\mathrm{M}_\odot\,\mathrm{yr}^{-1}) \f]

    where \f$Z\f$ is the metallicity, \f$\mathrm{SFE}\f$ is the star formation efficiency,
    \f$n_{\text{cl}}\f$ is the cloud number density, and \f$\dot{M}_*\f$ is the star formation rate.
    
    In SFRNormalized mode, this class assumes a constant star formation history
    over the past 10 or 30 Myr period (as determined by the sfrPeriod parameter). The model 
    properly accounts for the recollapse of gas clouds that occurs when stellar feedback is 
    insufficient to overcome gravity, resulting in multiple generations of star formation 
    within a single cloud. This recollapse contribution is pre-integrated over the time evolution
    (10 or 30 Myr) and cloud mass spectrum (\f$10^5\f$ to \f$10^{6.75}~\mathrm{M}_\odot\f$) with
    a power-law distribution (\f$dN/dM \propto M^{-1.8}\f$). T
    
    The model SEDs cover wavelengths from 0.01 to 3005 \f$\mu\mathrm{m}\f$ (UV through millimeter) and include stellar, 
    nebular, and dust continuum emission along with numerous emission lines from H II, PDR, and molecular 
    gas phases from Cloudy spectral synthesis calculations.
*/
class ToddlersSEDFamily : public SEDFamily
{
    /** The enumeration type indicating the SED calculation mode */
    ENUM_DEF(SedMode, Cloud, SFRNormalized)
        ENUM_VAL(SedMode, Cloud, "Individual cloud SEDs with time evolution")
        ENUM_VAL(SedMode, SFRNormalized, "SEDs normalized by star formation rate")
    ENUM_END()

    /** The enumeration type indicating the stellar template to use */
    ENUM_DEF(StellarTemplate, SB99Kroupa100Sin, BPASSChab100Bin, BPASSChab300Bin)
        ENUM_VAL(StellarTemplate, SB99Kroupa100Sin, "Starburst99 with Kroupa IMF (0.1-100 Msun) and single star evolution")
        ENUM_VAL(StellarTemplate, BPASSChab100Bin, "BPASS with Chabrier IMF (0.1-100 Msun) and binary star evolution")
        ENUM_VAL(StellarTemplate, BPASSChab300Bin, "BPASS with Chabrier IMF (0.1-300 Msun) and binary star evolution")
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

    ITEM_CONCRETE(ToddlersSEDFamily, SEDFamily,
                 "a TODDLERS SED family for emission from star-forming regions")
        PROPERTY_ENUM(sedMode, SedMode, "SED mode (Cloud or SFRNormalized)")
        ATTRIBUTE_DEFAULT_VALUE(sedMode, "SFRNormalized")

        PROPERTY_ENUM(stellarTemplate, StellarTemplate, "the stellar template, IMF, and evolution model to use")
        ATTRIBUTE_DEFAULT_VALUE(stellarTemplate, "SB99Kroupa100Sin")

        PROPERTY_BOOL(includeDust, "include dust processing in the SED models")
        ATTRIBUTE_DEFAULT_VALUE(includeDust, "true")

        PROPERTY_ENUM(resolution, Resolution, "the wavelength resolution")
        ATTRIBUTE_DEFAULT_VALUE(resolution, "Low")

        PROPERTY_ENUM(sfrPeriod, SFRPeriod, "the SFR integration time period (only for SFRNormalized mode)")
        ATTRIBUTE_DEFAULT_VALUE(sfrPeriod, "Period10Myr")
        ATTRIBUTE_RELEVANT_IF(sfrPeriod, "SFRNormalized")
    ITEM_END()

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family. The newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy, and its setup() function has been called. */
    explicit ToddlersSEDFamily(SimulationItem* parent, SedMode sedMode, bool includeDust, 
                               Resolution resolution, SFRPeriod sfrPeriod = SFRPeriod::Period10Myr);

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
    // Only one of these tables will be used, depending on the sedMode
    StoredTable<6> _cloudTable;       // 6D table: lambda, time, Z, SFE, n_cl, M_cl
    StoredTable<4> _sfrNormalizedTable; // 4D table: lambda, Z, SFE, n_cl
};

//////////////////////////////////////////////////////////////////////

#endif
