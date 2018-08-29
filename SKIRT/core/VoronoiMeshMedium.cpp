/*//////////////////////////////////////////////////////////////////
////     The SKIRT project -- advanced radiative transfer       ////
////       © Astronomical Observatory, Ghent University         ////
///////////////////////////////////////////////////////////////// */

#include "VoronoiMeshMedium.hpp"
#include "VoronoiMeshSnapshot.hpp"

////////////////////////////////////////////////////////////////////

Snapshot* VoronoiMeshMedium::createAndOpenSnapshot()
{
    // create and open the snapshot
    _voronoiMeshSnapshot = new VoronoiMeshSnapshot;
    _voronoiMeshSnapshot->open(this, filename(), "Voronoi sites");

    // configure the mass or density column (position columns are configured by the snapshot itself)
    switch (massType())
    {
    case MassType::MassDensity: _voronoiMeshSnapshot->importMassDensity(); break;
    case MassType::Mass: _voronoiMeshSnapshot->importMass(); break;
    case MassType::NumberDensity: _voronoiMeshSnapshot->importNumberDensity(); break;
    case MassType::Number: _voronoiMeshSnapshot->importNumber(); break;
    }

    // set the domain extent
    _voronoiMeshSnapshot->setExtent(domain());
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////

VoronoiMeshSnapshot* VoronoiMeshMedium::voronoiMesh() const
{
    return _voronoiMeshSnapshot;
}

////////////////////////////////////////////////////////////////////