# Neutral atmospheric-boundary-layer validation

This example provides the cases and diagnostics for validating the neutral ABL
implementation from issue #317 against the physical checks in issue #318. It
reproduces the uniform
`64 x 65 x 64`, aspect-ratio `pi` case from Deskos (2019) and the legacy
Incompact3d `input_neutral.i3d` configuration:

- `Lx x Ly x Lz = 3140 x 1000 x 3140 m`
- `z0 = 0.1 m`, `u_star = 0.45 m/s`
- pressure-gradient forcing `u_star^2 / delta`
- static Smagorinsky with Mason--Thomson wall damping `(Cs, n) = (0.14, 3)`
- periodic x/z, free-slip top, and rough-wall stress at the bottom

The matched legacy run uses `dt = 0.8`, 100000 steps, and accumulates its
reported statistics after step 25000. `input_pressure_gradient.x3d` uses the
same values. The companion `input_ekman.x3d` switches to geostrophic forcing
and Coriolis rotation to check the expected wind-direction veer.

## Running

Run each configuration in its own directory because the solver writes output
relative to its working directory:

```bash
mkdir -p runs/abl_pressure_gradient runs/abl_ekman

(cd runs/abl_pressure_gradient && \
  mpirun -np 1 ../../build/bin/xcompact \
  ../../examples/abl_validation/input_pressure_gradient.x3d)

(cd runs/abl_ekman && \
  mpirun -np 1 ../../build/bin/xcompact \
  ../../examples/abl_validation/input_ekman.x3d)
```

At every `n_output` interval after `profile_start_iter`, the ABL diagnostics
replace a small CSV with the latest horizontal/time average. The file includes
the mean velocity components, analytical log law, mean wall-stress vector,
and the friction velocity diagnosed from that stress. The pressure-gradient
run writes `abl_profile.csv`; the Ekman run writes `abl_ekman_profile.csv`.
Diagnostic accumulators are not stored in flow checkpoints: after a restart,
the CSV contains the average of samples collected by that restarted run.

Check the profiles and create `metrics.csv` and, when Matplotlib is available,
`profiles.png` with:

```bash
python3 examples/abl_validation/analyse.py \
  --profile runs/abl_pressure_gradient/abl_profile.csv \
  --ekman-profile runs/abl_ekman/abl_ekman_profile.csv \
  --check
```

The automated criteria require the surface-layer mean profile to remain close
to the analytical log law, diagnosed friction velocity within 5% of the
imposed value, outer-layer mean `Phi` within 20% of one, the expected
surface-layer overshoot within the published 1.2--1.5 band, and at least one
degree of Ekman veer. These thresholds are intentionally explicit in the
analysis command and can be overridden for sensitivity studies.

The raw checkpoints and flow fields are run artefacts and should not be
committed. Commit only the inputs, analysis outputs needed for review, and a
concise validation note.
