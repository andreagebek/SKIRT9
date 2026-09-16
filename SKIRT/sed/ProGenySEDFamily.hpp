/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROGENYSEDFAMILY_HPP
#define PROGENYSEDFAMILY_HPP

#include "SEDFamily.hpp"
#include "StoredTable.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the ProGenySEDFamily class represents a ProGeny family of single stellar population
    (SSP) %SED templates (Robotham and Bellstedt, 2024, RASTI, 4). We use the PG_Ch_Mi_Kn_Alpha_X version
    of the models with variable alpha-enhancement and a Chabrier IMF (0.1-100 \f$\mathrm{M}_\odot\f$).

    The SED templates are parametrized on metallicity (1e-6 - 0.06), age (0.1 Myr - 20 Gyr), and
    alpha-enhancement (-0.2 - 0.6), and scale with the initial mass of the SSP.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ M_\mathrm{init}\,(\mathrm{M}_\odot) \quad Z\,(\mathrm{dimensionless}) \quad
    t\,(\mathrm{yr}) \quad [\alpha/\mathrm{Fe}]\,(\mathrm{dimensionless}) \f] */

class ProGenySEDFamily : public SEDFamily
{
    ITEM_CONCRETE(ProGenySEDFamily, SEDFamily, "a ProGeny SED family for simple stellar populations")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked programmatically by classes that use a hard-coded SED
        family (as opposed to selected through the ski file). Before the constructor returns, the
        newly created object is hooked up as a child to the specified parent in the simulation
        hierarchy (so it will automatically be deleted), and its setup() function has been called.
        */
    explicit ProGenySEDFamily(SimulationItem* parent);

protected:
    /** This function opens the appropriate resource file (in SKIRT stored table format). */
    void setupSelfBefore() override;

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. It retrieves this
        range from the underlying stored table. */
    Range intrinsicWavelengthRange() const override;

    /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
        unit of wavelength) for the %SED with the specified parameters at the specified wavelength,
        or zero if the wavelength is outside of the %SED's intrinsic wavelength range. The number
        and type of parameters must match the information returned by the parameterInfo() function;
        if not the behavior is undefined. */
    double specificLuminosity(double wavelength, const Array& parameters) const override;

    /** This function constructs both the normalized probability density function (pdf) and the
        corresponding normalized cumulative distribution function (cdf) for the %SED with the
        specified parameters over the specified wavelength range. The function returns the
        normalization factor. The number and type of parameters must match the information returned
        by the parameterInfo() function; if not the behavior is undefined. */
    double cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
               const Array& parameters) const override;

    //====================== Data members =====================

private:
    StoredTable<4> _table;
};

////////////////////////////////////////////////////////////////////

#endif
