# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build Commands

```bash
# Configure (testing disabled by default)
cmake -B build -DENABLE_TESTING=ON

# Configure with Python bindings
cmake -B build -DENABLE_PYTHON=ON

# Configure with interactive visualization (requires SFML 2.5)
cmake -B build -DENABLE_INTERACTIVE=ON

# Build
cmake --build build

# Run all tests
ctest --test-dir build

# Run a specific test by name (Catch2 tag or test name)
./build/test/nugeometry-testsuite "[shape]"
./build/test/nugeometry-testsuite "test name substring"
```

External dependencies are fetched automatically via CPM (a CMake package manager wrapper around FetchContent). No manual dependency installation is needed.

## Architecture

NuGeometry is a C++17 geometry engine for neutrino detector simulations. It models detector geometry via a volume hierarchy, performs ray/sphere tracing for visualization, and computes interaction probabilities using material cross-sections.

### Object Hierarchy

```
World
└── PhysicalVolume   (LogicalVolume + Transform3D)
    └── LogicalVolume  (Shape + Material)
        ├── Shape      (Box, Sphere, Cylinder, Trapezoid, Polyhedra, CombinedShape)
        └── Material   (elements + fractions + density)
            └── Element  (Z, A, isotope composition)
```

**World** (`include/geom/World.hh`) is the top-level container. It owns a tree of `PhysicalVolume`s and is the main entry point for ray tracing and geometry queries.

**LogicalVolume** / **PhysicalVolume** (`include/geom/Volume.hh`) separate the shape/material definition from its placement. A `PhysicalVolume` is a positioned instance of a `LogicalVolume`.

### Shape System

All shapes (`include/geom/Shape.hh`) inherit from `Shape` and implement:
- `SignedDistance(point)` — negative inside, positive outside (SDF used for sphere tracing)
- `Intersect(ray)` — returns time parameter `t` along ray
- `GetBoundingBox()` — axis-aligned bounding box
- `Volume()` — analytic volume

**ShapeFactory** uses a registration pattern so GDML parsing can construct shapes by name at runtime:
```cpp
// Register a new shape type
ShapeFactory::Register<MyShape>("MyShape");
// Create by name (used by parser)
auto shape = ShapeFactory::Initialize("Box", node);
```
`RegistrableShape<Derived>` is the CRTP base that handles registration.

`CombinedShape` implements boolean CSG operations (union, intersection, subtraction) on two child shapes.

### GDML Parser

`GDMLParser` (`include/geom/Parser.hh`, `src/geom/Parser.cc`) reads GDML XML files (the Geant4 geometry format) using pugixml. It constructs the full World hierarchy including materials, shapes, and volume placements.

### Material / Element System

`Material` stores element compositions by mass fraction or atom count, and computes number densities used for cross-section calculations. `Element` tracks Z, A, and optional isotope breakdowns. A static registry in `Element` provides common elements by symbol.

### Physics / Simulation

`DetectorSim` (`include/geom/DetectorSim.hh`) generates neutrino interaction events: it samples interaction points along a ray through detector volumes, weighting by material cross-sections and number densities.

`TestGen` / `ExProb` provide probability calculation utilities used by the `prob_test` and `prob_test2` executables.

### Ray Tracing / Visualization

`main.cc` renders detector geometry as a PPM image using sphere tracing (iterative SDF stepping) with Phong shading and shadows. The camera (`include/geom/Camera.hh`) projects rays through the scene.

### Python Bindings

When built with `-DENABLE_PYTHON=ON`, pybind11 compiles `PythonInterface.cc` into a shared module `nugeom` that exposes geometry creation and querying to Python.

## Key Conventions

- Headers use `.hh`, sources use `.cc`
- Code style is enforced by `.clang-format`
- Pre-commit hooks are configured in `.pre-commit-config.yaml`
- Tests use Catch2 v2 with Trompeloeil for mocking; all tests live in a single executable `nugeometry-testsuite`
- The `TESTING` preprocessor define is set when `ENABLE_TESTING=ON`
