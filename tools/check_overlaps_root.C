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
#include "TGeoMaterial.h"
#include "TGeoNavigator.h"
#include "TGeoNode.h"
#include "TGeoVolume.h"
#include "TObjArray.h"
#include "TSystem.h"
#include <cstdio>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

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

// Step ROOT's TGeoNavigator along each ray in tools/compare_rays.txt and write
// tools/root_segments.txt in the same "rayIndex  len:mat ..." format the
// NuGeometry [navcompare] test emits, so the two traversals can be diffed.
void navigate_rays(const char *gdml_path = "nd_hall_with_lar_tms_nosand.gdml",
                   const char *rays_file = "tools/compare_rays.txt",
                   const char *out_file = "tools/root_segments.txt") {
    gSystem->Load("libGeom");
    TGeoManager *geo = TGeoManager::Import(gdml_path);
    if(!geo) {
        printf("ERROR: could not import '%s'\n", gdml_path);
        return;
    }
    std::ifstream in(rays_file);
    if(!in) {
        printf("ERROR: could not open '%s'\n", rays_file);
        return;
    }
    std::ofstream out(out_file);
    out.precision(10);

    std::string line;
    int idx = 0;
    while(std::getline(in, line)) {
        if(line.empty()) continue;
        std::istringstream ss(line);
        double ox, oy, oz, dx, dy, dz;
        if(!(ss >> ox >> oy >> oz >> dx >> dy >> dz)) continue;

        TGeoNavigator *nav = geo->GetCurrentNavigator();
        if(!nav) nav = geo->AddNavigator();
        nav->InitTrack(ox, oy, oz, dx, dy, dz);

        out << idx;
        int guard = 0;
        while(!nav->IsOutside() && guard++ < 100000) {
            TGeoVolume *vol = nav->GetCurrentVolume();
            const char *mat =
                (vol && vol->GetMaterial()) ? vol->GetMaterial()->GetName() : "UNKNOWN";
            nav->FindNextBoundaryAndStep();
            const double step = nav->GetStep();
            out << "  " << step << ":" << mat;
        }
        out << "\n";
        ++idx;
    }
    printf("Wrote %d ray segmentations to %s\n", idx, out_file);
}
