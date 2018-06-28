/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef IMPORTEDSOURCE_HPP
#define IMPORTEDSOURCE_HPP

#include "Source.hpp"
#include "Array.hpp"
#include "SED.hpp"
class Snapshot;

//////////////////////////////////////////////////////////////////////

/** ImportedSource is an abstract class representing a primary radiation source with a spatial and
    spectral luminosity distribution imported from an input file. The input data is usually derived
    from a hydrodynamical simulation snapshot. Various types of snapshots are supported by
    subclasses of this class. Refer to the subclass documentation for information on the file
    format.

    Usually, the input file defines a spatial distribution through smoothed particles, which must
    be interpolated and summed, or through adjacent cells that partition the spatial domain. At the
    level of this abstract class, we use the generic term \em entity for referring to either a
    particle or a cell.

    In addition to spatial information, each entity in the snapshot carries properties that allow
    selecting a particular %SED from a parameterized %SED family. The present class requires the
    user to configure an SEDFamily object for this purpose. The number, type, and order of
    parameters is defined by the %SED family. For each entity, the %SED family is requested to
    select and properly scale a specific %SED based on the entity's properties. Combining the
    spatial and spectral information for an entity yields its contribution to the imported
    radiation source.

    The input file may also include a seperate (bulk) velocity vector for each entity. When this
    option is enabled, the appropriate Doppler shift is taken into account when launching photon
    packets. Apart from the anisotropy resulting from this optional Doppler shift, the radiation
    emitted by this primary source is always isotropic. It is also always unpolarized. */
class ImportedSource : public Source
{
    ITEM_ABSTRACT(ImportedSource, Source, "a primary source imported from snapshot data")

    PROPERTY_STRING(filename, "the name of the file to be imported")

    ATTRIBUTE_SUB_PROPERTIES_HERE()

    PROPERTY_BOOL(importVelocity, "import velocity components (3 columns)")
        ATTRIBUTE_DEFAULT_VALUE(importVelocity, "false")

    PROPERTY_ITEM(sed, SED, "the spectral energy distribution for the source XXX")
        ATTRIBUTE_DEFAULT_VALUE(sed, "SunSED")

    ITEM_END()

    //============= Construction - Setup - Destruction =============

protected:
    /** This function imports the snapshot data from the input file through a Snapshot object of
        the appropriate type. Specifically, it first calls the createSnapshot() function, which
        must be implemented in a subclass, to construct and open a Snapshot object of the
        appropriate type. It then passes the user-configurable options of this class to the
        Snapshot object and tells it to import the data.

        Finally, the function constructs a vector with the luminosities (integrated over the
        primary source wavelength range) for all imported entities. This information is used when
        deciding how many photon packets should be launched from each entity. */
    void setupSelfBefore() override;

    /** This function constructs a new Snapshot object of the type appropriate for the subclass,
        calls its open() function, and returns a pointer to the object. Ownership of the Snapshot
        object is transferred to the caller. */
    virtual Snapshot* createAndOpenSnapshot() = 0;

    /** The destructor deletes the snapshot object, if present. */
    ~ImportedSource();

    //======================== Other Functions =======================

public:
    /** This function returns the dimension of the source, which is always 3 for an imported
        source. */
    int dimension() const override;

    /** This function returns the luminosity \f$L\f$ (i.e. radiative power) of the source
        integrated over the wavelength range of primary sources (configured for the source system
        as a whole) and across its complete spatial domain. */
     double luminosity() const override;

     /** This function returns the specific luminosity \f$L_\lambda\f$ (i.e. radiative power per
         unit of wavelength) of the source at the specified wavelength, or zero if the wavelength is
         outside the wavelength range of primary sources (configured for the source system as a
         whole) or if the source simply does not emit at the wavelength. */
     double specificLuminosity(double wavelength) const override;

     /** This function performs some preparations for launching photon packets. It is called in
         serial mode before each segment of photon packet launches, providing the history indices
         mapped by the source system to this particular source. See the description of the
         SourceSystem class for more background information.

         This function distributes the provided range of history indices over the individual
         entities imported by this source, creating a map for use when actually launching the
         photon packets. The number of photon packets allocated to each entity is determined as
         follows:

         \f[ N_m = \left[ (1-\xi) \frac{L_m}{L} + \xi \frac{1}{M} \right] N_s \f]

         where \f$N_s\f$ is the total number of photon packets to be launched by this source,
         \f$N_m\f$ is the number of photon packets to be launched by entity \f$m\f$, \f$L_m\f$ is
         the luminosity of source \f$m\f$, \f$L\f$ is the total luminosity for this source, \f$M\f$
         is the number of entities in this source, and \f$\xi\f$ is the \em emissionBias property
         value of the source system. */
     void prepareForLaunch(double sourceBias, size_t firstIndex, size_t numIndices) override;

     /** This function causes the photon packet \em pp to be launched from the source using the
         given history index and luminosity contribution. XXX. */
    void launch(PhotonPacket* pp, size_t historyIndex, double L) const override;

    //======================== Data Members ========================

private:
    // initialized during setup
    Snapshot* _snapshot{nullptr};
    double _L{0};       // the total bolometric luminosity of all entities (absolute number)
    Array _Lv;          // the relative bolometric luminosity of each entity (normalized to unity)

    // intialized by prepareForLaunch()
    Array _Wv;          // the relative launch weight for each entity (normalized to unity)
    vector<size_t> _Iv; // first history index allocated to each entity (with extra entry at the end)
};

//////////////////////////////////////////////////////////////////////

#endif