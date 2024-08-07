/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef SPINFLIPSEDFAMILY_HPP
#define SPINFLIPSEDFAMILY_HPP

#include "SEDFamily.hpp"

//////////////////////////////////////////////////////////////////////

/** An instance of the SpinFlipSEDFamily class represents a family of Gaussian spectra around the
    central 21-cm spin-flip wavelength \f$\lambda_\mathrm{sf}\f$, reflecting the thermal sub-grid
    motion in the source. The %SED family is parameterized on the bolometric 21-cm line luminosity
    \f$L_\mathrm{sf}\f$ and the spectral dispersion \f$s\f$ in velocity units. Using the photon
    velocity shift \f[ v = \frac{\lambda - \lambda_\mathrm{sf}} {\lambda_\mathrm{sf}} \,c \f] as
    the spectral variable, the spectrum can be written as (see also the LyaGaussianSED class): \f[
    L_v(v) = \frac{L_\mathrm{sf}}{s\,\sqrt{2\pi}} \,\exp\left( -\frac{v^2}{2s^2} \right). \f]

    The intrinsic range for the complete %SED family is taken to be approximately \f$\pm 9s\f$
    around the center for a dispersion of \f$s=1000\,\mathrm{km/s}\f$. This results in a range of
    approximately \f$20.47 \mathrm{cm} \le \lambda \le 21.74 \mathrm{cm}\f$. The source wavelength
    range configured by the user must fully contain this intrinsic wavelength range.

    Whenever a tabular form of a Gaussian %SED is requested, this class uses 100 wavelength points
    per dispersion unit on a regular linear grid.

    When imported from a text column file, the parameters for this %SED family must appear in the
    following order in the specified default units (unless these units are overridden by column
    header info): \f[ L_\mathrm{sf}\,(\mathrm{W}) \quad s\,(\mathrm{km/s}) \f] */
class SpinFlipSEDFamily : public SEDFamily
{
    ITEM_CONCRETE(SpinFlipSEDFamily, SEDFamily, "a family of Gaussian spectra around the central spin-flip wavelength")
    ITEM_END()

    //====================== Other functions =====================

public:
    /** This function returns the number and type of parameters used by this particular %SED family
        as a list of SnapshotParameter objects. Each of these objects specifies unit information
        and a human-readable descripton for the parameter. */
    vector<SnapshotParameter> parameterInfo() const override;

    /** This function returns the intrinsic wavelength range of the %SED family. For the
        SpinFlipSEDFamily, the intrinsic range is determined as decribed in the class header. */
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
};

////////////////////////////////////////////////////////////////////

#endif