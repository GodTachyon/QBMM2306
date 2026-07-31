# PbeBubbleFoam: Use Sauter Mean Diameter (d32) for Interfacial Forces

## Objective

Modify the `pbeBubbleFoam` solver (OpenQBMM, OpenFOAM v2306) so that the
interfacial forces are calculated using the **Sauter mean diameter** `d32`
instead of the per-node (per-size-class) bubble diameters.

## Summary of behaviour after the change

- **Interfacial forces** — drag, lift, wall lubrication, turbulent dispersion,
  and the `Re`/`We`/`Eo`/`EoH`/`Ta` correlations that feed them (plus wall
  damping and aspect-ratio models) — are evaluated with the Sauter mean
  diameter `d32`.
- **Momentum construction** — the B&G bubble-pressure viscosity in
  `divDevRhoReff1/2` keeps using the alpha-weighted mean diameter `d_`
  (unchanged behaviour). `d()` was deliberately left **non-virtual** so that
  d32 does not leak into the momentum equations.
- `d32_` is computed every time step in `polydispersePhaseModel::correct()` as
  the ratio of the third and second size moments and is written out as a field.

## Files changed

Only one source file differs from the repository baseline:

`OpenQBMM-OpenFOAM-v2306/applications/solvers/multiphase/twoPhaseSystemPbeBubble/phaseModels/polydispersePhaseModel/polydispersePhaseModel.H`

### 1. Override `ds(label nodei)` to return `d32_`

All blended interfacial force models obtain the bubble diameter through the
virtual accessor `pair_.dispersed().ds(nodei)`. Overriding it to return the
Sauter mean diameter for every node routes every diameter-based force model
through `d32_` without touching the shared `interfacialModels` library:

```cpp
//- Return the diameter for nodei
//  Override to return the Sauter mean diameter so that all
//  interfacial force models evaluate using d32
virtual const volScalarField& ds(const label nodei) const
{
    return d32_;
}
```

Previously this returned the per-node diameter `ds_[nodei]`.

### 2. Remove the now-inert `d()` override

The base-class `phaseModel::d()` is **non-virtual** and is only used by the
momentum construction (B&G viscosity in `divDevRhoReff1/2`). A `d()` override
returning `d32_` in the derived class is dead code (never dispatched through a
`const phaseModel&`) and would be misleading, so it was removed:

```cpp
// Removed:
// virtual const volScalarField& d() const { return d32_; }
```

> Note: `twoPhaseSystemPbeBubble/phaseModels/phaseModel/phaseModel.H` was
> temporarily edited to make `d()` virtual and later reverted. It now matches
> the repository baseline exactly (non-virtual `d()` returning `d_`).

## Why this works

- Every drag/lift/wall-lubrication/turbulent-dispersion model and every
  dimensionless correlation in `phasePair` queries the diameter via
  `pair_.dispersed().ds(nodei)` (virtual) or `pair_.dispersed().d()`
  (non-virtual).
- The virtual `ds(nodei)` override dispatches to `d32_` for all interfacial
  force evaluations.
- Drag collapses to the monodisperse-with-d32 form:
  `Kd = Σᵢ alphas(i)·Ki(d32) ≈ α·Ki(d32)` with
  `Ki = 0.75·CdRe·Cs·ρc·νc/d32²`.
- The non-virtual `d()` keeps the B&G bubble-pressure viscosity
  (`twoPhaseSystemPbeBubble.C`, `divDevRhoReff1/2`) on the alpha-weighted mean
  diameter `d_`, leaving momentum construction/solving unaffected.
- Per-node `ds_` fields are still maintained internally (used to build
  `alphas_` and to accumulate `d32_`); coalescence/breakup kernels use
  quadrature-node diameters directly and are unaffected.

## Git diff

```diff
diff --git a/OpenQBMM-OpenFOAM-v2306/applications/solvers/multiphase/twoPhaseSystemPbeBubble/phaseModels/polydispersePhaseModel/polydispersePhaseModel.H b/OpenQBMM-OpenFOAM-v2306/applications/solvers/multiphase/twoPhaseSystemPbeBubble/phaseModels/polydispersePhaseModel/polydispersePhaseModel.H
index 0e52546..7ad75ec 100644
--- a/OpenQBMM-OpenFOAM-v2306/applications/solvers/multiphase/twoPhaseSystemPbeBubble/phaseModels/polydispersePhaseModel/polydispersePhaseModel.H
+++ b/OpenQBMM-OpenFOAM-v2306/applications/solvers/multiphase/twoPhaseSystemPbeBubble/phaseModels/polydispersePhaseModel/polydispersePhaseModel.H
@@ -207,9 +207,11 @@ public:
         }
 
         //- Return the diameter for nodei
+        //  Override to return the Sauter mean diameter so that all
+        //  interfacial force models evaluate using d32
         virtual const volScalarField& ds(const label nodei) const
         {
-            return ds_[nodei];
+            return d32_;
         }
 	/* [Removed since this is an override from phaseModel.H, we need the default for now]
         //- Return the velocity for nodei
@@ -236,9 +238,6 @@ public:
         //- Correct the phase diameter to mean diameter from moments [NEW ADDITION]
         //virtual void correctDiameter();
         
-        // Override base class accessor, for computing the sauter mean diameter [NEW ADDITION]
-	    virtual const volScalarField& d() const { return d32_; }
-
         //- Correct the phase properties
         virtual void correct();
```

## Build

Rebuild only the dynamic library via the local wmake (the `pbeBubbleFoam`
executable is dynamically linked against it):

```sh
cd applications/solvers/multiphase
wmake twoPhaseSystemPbeBubble        # builds libpdPhaseSystemPbe.so
```

Build result: **no errors, no warnings**.
Output library:
`/students/2024/ashok/OpenFOAM/ashok-v2306/platforms/linux64GccDPInt32Opt/lib/libpdPhaseSystemPbe.so`

The executable does not need a rebuild: its only direct virtual calls on
`phaseModel` (`nNodes()`, `alphas()`) are declared before `d()` and are
therefore unaffected by the `d()` virtual/non-virtual layout.

## Verification

- Interfacial force coefficients (`Kd`, `F`, `D`) are now based on `d32` (all
  quadrature nodes share the same Sauter mean diameter).
- `divDevRhoReff1/2` (B&G viscosity) still uses the alpha-weighted mean
  diameter `d_`.
- Run the case and check the `d32` field output and the `Kd`/drag behaviour.
