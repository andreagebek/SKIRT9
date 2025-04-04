/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "SpinFlipSEDFamily.hpp"
#include "Constants.hpp"
#include "NR.hpp"

//////////////////////////////////////////////////////////////////////

namespace
{
    // number of wavelength points in tabulated result per dispersion unit
    const int numWavelengthsPerDispersionUnit = 100;

    // Gaussian centered on 0 with dispersion of 1, evaluated at x
    double unitGaussian(double x)
    {
        return (0.5 * M_SQRT1_2 * M_2_SQRTPI) * exp(-0.5 * x * x);
    }
}

////////////////////////////////////////////////////////////////////

vector<SnapshotParameter> SpinFlipSEDFamily::parameterInfo() const
{
    return {SnapshotParameter::custom("line luminosity", "bolluminosity", "W"),
            SnapshotParameter::custom("dispersion", "velocity", "km/s")};
}

////////////////////////////////////////////////////////////////////

Range SpinFlipSEDFamily::intrinsicWavelengthRange() const
{
    constexpr Range range(Constants::lambdaSpinFlip() * (1. - 0.03), Constants::lambdaSpinFlip() * (1. + 0.03));
    return range;
}

////////////////////////////////////////////////////////////////////

double SpinFlipSEDFamily::specificLuminosity(double wavelength, const Array& parameters) const
{
    double L = parameters[0];
    double s = parameters[1];
    double wavelengthCenter = Constants::lambdaSpinFlip();
    double wavelengthDispersion = s * wavelengthCenter / Constants::c();
    return L * unitGaussian((wavelength - wavelengthCenter) / wavelengthDispersion) / wavelengthDispersion;
}

////////////////////////////////////////////////////////////////////

double SpinFlipSEDFamily::cdf(Array& lambdav, Array& pv, Array& Pv, const Range& wavelengthRange,
                              const Array& parameters) const
{
    double L = parameters[0];
    double s = parameters[1];
    double wavelengthCenter = Constants::lambdaSpinFlip();
    double wavelengthDispersion = s * wavelengthCenter / Constants::c();

    // build an appropriate grid
    size_t n = numWavelengthsPerDispersionUnit * wavelengthRange.width() / wavelengthDispersion;
    NR::buildLinearGrid(lambdav, wavelengthRange.min(), wavelengthRange.max(), n);

    // calculate the tabulated values
    pv.resize(n + 1);
    for (size_t i = 0; i <= n; ++i)
        pv[i] = L * unitGaussian((lambdav[i] - wavelengthCenter) / wavelengthDispersion) / wavelengthDispersion;

    // calculate the cumulative distribution and normalization
    return NR::cdf2(false, lambdav, pv, Pv);
}

////////////////////////////////////////////////////////////////////