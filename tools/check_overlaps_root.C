// External geometry overlap cross-check using ROOT's TGeoManager.
//
// Independently imports the GDML and runs ROOT's own overlap checker.  If ROOT
// flags the same volume pairs that NuGeometry's sweep/containment scan flags
// (World::CheckSweepConsistency / the [overlapscan] test), the overlaps are
// real GDML errors, not a bug in NuGeometry's parser/intersection code.
//
// Usage:
//   root -l -b -q 'tools/check_overlaps_root.C("nd_hall_with_lar_tms_nosand.gdml")'
//
// Notes:
//   * precision 1e-4 cm (1 micron) matches the sweep's prune threshold.
//   * runs both the mesh/bounding-box check and the "s" sampling check
//     (the latter catches overlaps the surface mesh misses).

#include "TGeoManager.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <cstdio>

void check_overlaps_root(const char *gdml_path = "nd_hall_with_lar_tms_nosand.gdml") {
    gSystem->Load("libGeom"); // ROOT geometry (TGeoManager) + GDML import
    TGeoManager *geo = TGeoManager::Import(gdml_path);
    if(!geo) {
        printf("ERROR: could not import '%s'\n", gdml_path);
        return;
    }
    const double prec = 1e-4; // cm
    printf("Imported '%s': %d logical volumes, %d total nodes in tree\n", gdml_path,
           geo->GetListOfVolumes() ? geo->GetListOfVolumes()->GetEntries() : -1, geo->CountNodes());
    // Spot-check the deep nesting our scan flags (volTPCActive should sit inside
    // the cryostat/detector chain, not overlap volHalfDetector).
    for(const char *vn : {"volTPCActive", "volHalfDetector", "volWLS", "volOpticalDet"}) {
        TGeoVolume *v = geo->GetVolume(vn);
        printf("  volume '%s': %s, %d daughters\n", vn, v ? "present" : "MISSING",
               v ? v->GetNdaughters() : -1);
    }

    printf("\n==== Mesh/bounding-box overlap check (precision %g cm) ====\n", prec);
    geo->CheckOverlaps(prec);
    geo->PrintOverlaps();

    printf("\n==== Sampling overlap check (precision %g cm, option 's') ====\n", prec);
    geo->CheckOverlaps(prec, "s");
    geo->PrintOverlaps();

    TObjArray *ov = geo->GetListOfOverlaps();
    printf("\nTotal overlaps reported by ROOT: %d\n", ov ? ov->GetEntriesFast() : 0);
}
