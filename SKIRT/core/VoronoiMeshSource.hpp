/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef VORONOIMESHSOURCE_HPP
#define VORONOIMESHSOURCE_HPP

#include "MeshSource.hpp"

////////////////////////////////////////////////////////////////////

/** A VoronoiMeshSource instance represents a primary radiation source with a spatial and spectral
    luminosity distribution described by a list of sites generating a Voronoi tesselation of a
    cubodail domain. The data is usually extracted from a cosmological simulation snapshot, and it
    must be provided in a column text file formatted as described below.

    Refer to the description of the TextInFile class for information on overall formatting and on
    how to include header lines specifying the units for each column in the input file. In case the
    input file has no unit specifications, the default units mentioned below are used instead. The
    number of columns expected in the input file depends on the options configured by the user for
    this VoronoiMeshSource instance, including the selected %SEDFamily.

    \f[ x\,(\mathrm{pc}) \quad y\,(\mathrm{pc}) \quad z\,(\mathrm{pc}) \quad [ v_x\,(\mathrm{km/s})
    \quad v_y\,(\mathrm{km/s}) \quad v_z\,(\mathrm{km/s}) \quad [ \sigma_v\,(\mathrm{km/s}) ] ]
    \quad [M_\mathrm{curr}\,(\mathrm{M}_\odot)] \quad b\,(1) \quad \dots \text{SED family parameters}
    \dots \f]

    The first three columns are the \f$x\f$, \f$y\f$ and \f$z\f$ coordinates of the Voronoi site
    (i.e. the location defining a particular Voronoi cell). If the \em importVelocity option is
    enabled, the next three columns specify the \f$v_x\f$, \f$v_y\f$, \f$v_z\f$ bulk velocity
    components of the source population represented by the cell corresponding to the site. If
    additionally the \em importVelocityDispersion option is enabled, the next column specifies the
    velocity dispersion \f$\sigma_v\f$, adjusting the velocity for each photon packet launch with a
    random offset sampled from a spherically symmetric Gaussian distribution. If the
    \em importCurrentMass option is enabled, the next column provides the current mass of the
    particle, \f$M_\mathrm{curr}\f$. This mass is currently only used for probing the input model.
    If the \em importBias option is enabled, the next column specifies the bias parameter, \f$b\f$,
    which is used to bias the photon sampling for each cell (see the documentation of the
    ImportedSource class).

    The remaining columns specify the parameters required by the configured %SED family to select
    and scale the appropriate %SED. For example for the Bruzual-Charlot %SED family, the remaining
    columns provide the initial mass, the metallicity, and the age of the stellar population
    represented by the cell corresponding to the site. Refer to the documentation of the configured
    type of SEDFamily for information about the expected parameters and their default units. */
class VoronoiMeshSource : public MeshSource
{
    ITEM_CONCRETE(VoronoiMeshSource, MeshSource, "a primary source imported from data represented on a Voronoi mesh")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function constructs a new VoronoiMeshSnapshot object, calls its open() function,
        passes it the domain extent configured by the user (using properties offered by the
        MeshSource base class), and returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    Snapshot* createAndOpenSnapshot() override;
};

////////////////////////////////////////////////////////////////////

#endif
