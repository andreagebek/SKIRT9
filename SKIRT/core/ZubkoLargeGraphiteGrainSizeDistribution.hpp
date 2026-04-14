/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef ZUBKOLARGEGRAPHITEGRAINSIZEDISTRIBUTION_HPP
#define ZUBKOLARGEGRAPHITEGRAINSIZEDISTRIBUTION_HPP

#include "GrainSizeDistribution.hpp"

////////////////////////////////////////////////////////////////////

/** ZubkoLargeGraphiteGrainSizeDistribution represents the dust grain size distribution and grain size
    range for the large graphite population in model BARE_GR_S of Zubko, Dwek & Arendt (2004, ApJS, 152,
    211). The size distribution function is scaled to obtain the appropriate dust mass per hydrogen
    atom for the graphite grain population. */
class ZubkoLargeGraphiteGrainSizeDistribution : public GrainSizeDistribution
{
    ITEM_CONCRETE(ZubkoLargeGraphiteGrainSizeDistribution, GrainSizeDistribution,
                  "a Zubko, Dwek & Arendt size distribution for large graphite dust grains")
    ITEM_END()

    //============= Construction - Setup - Destruction =============

public:
    /** This constructor can be invoked by classes that wish to hard-code the creation of a new
        grain size distribution object of this type (as opposed to creation through the ski file).
        Before the constructor returns, the newly created object is hooked up as a child to the
        specified parent in the simulation hierarchy (so it will automatically be deleted), and its
        setup() function has been called. */
    explicit ZubkoLargeGraphiteGrainSizeDistribution(SimulationItem* parent);

    //======================== Other Functions =======================

public:
    /** This function returns the hard-coded minimum grain size \f$a_\text{min}\f$ for this grain
        size distribution. */
    double amin() const override;

    /** This function returns the hard-coded maximum grain size \f$a_\text{max}\f$ for this grain
        size distribution. */
    double amax() const override;

    /** This function returns the hard-coded value of this grain size distribution. */
    double dnda(double a) const override;
};

////////////////////////////////////////////////////////////////////

#endif
