/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef GRAINSIZEDISTRIBUTION_HPP
#define GRAINSIZEDISTRIBUTION_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** GrainSizeDistribution is an abstract base class that represents a size distribution for the
    dust grains in a particular dust population. Specifically, it represents a function
    \f$\text{dnda}(a)\f$ that specifies the relative number of dust grains with size \f$a\f$ in the
    population, \f[ \text{dnda}(a) \propto \frac{\text{d}n_\text{D}}{\text{d}a} \qquad
    \text{for}\quad a_\text{min} \leq a \leq a_\text{max}. \f] The function is scaled arbitrarily;
    an appropriate proportionality factor is determined elsewhere by specifying a normalization for
    the amount of dust in the population.

    The GrainSizeDistribution class offers access to the size distribution range and the size
    distribution value within that range. It expects each subclass to implement the functions
    declared in this interface, i.e. the functions amin() and amax() to specify the grain size
    range, and the function dnda() to specify the grain size distribution function within that
    range. */
class GrainSizeDistribution: public SimulationItem
{
    ITEM_ABSTRACT(GrainSizeDistribution, SimulationItem, "a dust grain size distribution")
    ITEM_END()

public:
    /** This function returns the minimum grain size \f$a_\text{min}\f$, i.e. the lower limit of
        the distribution. */
    virtual double amin() const = 0;

    /** This function returns the maximum grain size \f$a_\text{max}\f$, i.e. the upper limit of
        the distribution. */
    virtual double amax() const = 0;

    /** This function returns the value of the distribution \f$\text{dnda} \propto
        \frac{\text{d}n_\text{D}}{\text{d}a}\f$ for a given grain size \f$a\f$. If
        \f$a<a_\text{min}\f$ or \f$a>a_\text{max}\f$ the result is undefined (i.e. the function
        does not need to check the bounds). */
    virtual double dnda(double a) const = 0;
};

////////////////////////////////////////////////////////////////////

#endif
