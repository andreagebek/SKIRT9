/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#ifndef PROBE_HPP
#define PROBE_HPP

#include "SimulationItem.hpp"

////////////////////////////////////////////////////////////////////

/** Probe is an abstract class representing probes that output information on internal simulation
    data before or during the simulation run. Refer to the ProbeSystem class for more information.

    Probe subclasses \em must adhere to the following rules. In this discussion, \em target refers
    to the item(s) in the simulation hierarchy from which the probe retrieves information.
    - Data encapsulation: a probe can use public interfaces only, even if these interfaces are in
      some cases written specifically to support the probe (refer to the scenarios described in
      the ProbeSystem class).
    - Read-only access: a probe cannot (cause to) change the data structures held by the target.
      There are two exceptions to this rule: (1) during its own setup, a probe can cause setup of
      a target through the find() or interface() functions, and (2) a probe can use a public
      function, provided by the target for this purpose, to install a call-back function that will
      be invoked by the target.
    - Interprobe independence: a probe cannot look for or depend on another probe, nor on the order
      of the various probes in the list held by the probe system.
*/
class Probe : public SimulationItem
{
    ITEM_ABSTRACT(Probe, SimulationItem, "a probe")
        PROPERTY_STRING(probeName, "the name for this probe")
    ITEM_END()

    //============== Functions to be implemented in subclass =============

protected:
    /** This enumeration indicates when to perform probing: after setup or after the complete
        simulation run. */
    enum class When { Setup, Run };

    /** This function returns an enumeration indicating when probing should be performed for this
        probe. The default implementation in this base class returns \c Setup. A subclass needs to
        override this function only if it (may) require probing at a different time. */
    virtual When when() const;

    /** This function is called after the simulation has been fully setup but before the probe()
        function is called. It can implemented by a subclass that needs to perform some
        initialization that requires the simulation tom be fully setup. The default implementation
        in this base class does nothing. */
    virtual void initialize();

    /** This function must be implemented in each subclass to produce the relevant probing output.
        It is called either at the end of setup or at the end of the simulation run, depending on
        the return value of the when() function. */
    virtual void probe() = 0;

    //=================== Functions implemented here ===================

public:
    /** This function returns the probe name as human-readable name for the simulation item, so
        that it can be used in log messages to identify the probe and differentiate it from other
        probes. */
    string itemName() const override;

    /** This function is called at the end of the setup phase, i.e. after all simulation items have
        performed setup. It first calls initialize(), and if when() returns \c Setup, it then calls
        probe(). */
    void probeSetup();

    /** This function is called at the end of the run phase, i.e. after all photon packets have
        been emitted and detected. If when() returns \c Run, this function calls probe(); otherwise
        it does nothing. */
    void probeRun();
};

////////////////////////////////////////////////////////////////////

#endif
