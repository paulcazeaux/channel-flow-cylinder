# Benchmark Targets

## Reference: Schafer & Turek (1996) DFG Benchmark

Schafer, M. & Turek, S. (1996). "Benchmark computations of laminar flow around a cylinder."
In: Flow Simulation with High-Performance Computers II, Notes on Numerical Fluid Mechanics 52,
Vieweg, pp. 547–566.

## Domain & Parameters

```
Channel:     [0, 2.2] × [0, 0.41] m
Cylinder:    center (0.2, 0.2), radius r = 0.05 m, diameter D = 0.1 m
Viscosity:   ν = 1e-3 m²/s
Density:     ρ = 1.0 kg/m³
Inlet:       u(0,y) = 4·U_max·y·(H-y)/H², U_max = 0.3 m/s
Mean velocity: U_mean = (2/3)·U_max = 0.2 m/s
Re = U_mean·D/ν = 0.2·0.1/1e-3 = 20   (Schafer-Turek 2D-1)
```

## Schafer-Turek 2D-1 Benchmark Values (Re = 20, steady)

| Quantity | Reference interval | Notes |
|----------|--------------------|-------|
| C_D (drag coefficient) | 5.5700 – 5.5900 | Primary validation target |
| C_L (lift coefficient) | 0.0104 – 0.0110 | Near-zero (symmetric eddies) |
| Δp (pressure drop, front–rear of cylinder) | 0.1172 – 0.1176 Pa | Secondary check |

Coefficients defined as:
  C_D = 2·F_D / (ρ·U_mean²·D·1)
  C_L = 2·F_L / (ρ·U_mean²·D·1)

## This Project: Re = 5

At Re=5 the flow is more viscous; quantitative reference values are not tabulated in
Schafer-Turek (which targets Re=20 and Re=100). Expected behavior:
- C_D > 5.57 (drag increases as Re decreases in Stokes-like regime; roughly C_D ∝ 1/Re)
- C_L ≈ 0 (flow remains symmetric and steady at Re=5)
- Closed symmetric recirculation eddies behind cylinder
- No vortex shedding (onset ≈ Re=47 for cylinder)

To run the Re=20 case for validation: set U_max = 0.3 m/s in `run_simulation.py`
and compare C_D against the interval above.

## Validation Strategy

1. Run Re=20 case → check C_D ∈ [5.57, 5.59]
2. Run Re=5 case → verify C_L ≈ 0, symmetric velocity field
3. Inspect velocity profiles at x=0.5 m cross-section for expected parabolic recovery
