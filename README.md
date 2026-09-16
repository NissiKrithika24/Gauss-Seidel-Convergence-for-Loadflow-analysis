# Gauss–Seidel Load Flow Solver (MATLAB)

A general-purpose Gauss–Seidel power-flow (load-flow) solver. It builds the
Y-bus admittance matrix from user-supplied line data, accepts Slack/PV/PQ
bus data, and iteratively solves for node voltages until convergence, with
optional successive over-relaxation and PV-bus reactive-power (Q) limit
handling.

## Files

| File | Description |
|---|---|
| `gauss_seidel_loadflow.m` | Original interactive script — prompts for bus/line data via `input()`. Run this to solve your own system. |
| `gauss_seidel_loadflow_demo.m` | Non-interactive version with a hardcoded 3-bus sample system, so it can run unattended (used to generate the results below). |
| `sample_input.txt` | The exact values to type in when running the interactive script, for the 3-bus sample system. |
| `gauss_seidel_output.txt` | Full console output from running the demo script. |
| `bus_results.csv` | Final bus voltage magnitudes and angles. |
| `convergence_history.csv` | Max voltage change per iteration (for plotting convergence). |
| `convergence_plot.png` | Semilog plot of convergence history. |
| `gauss_seidel_results.mat` | Full MATLAB workspace (Y-bus, voltages, bus data, etc.) for reuse. |

## How it works

1. **Y-bus formation** — off-diagonal `Y(i,j) = -1/(R+jX)` for each line;
   diagonal `Y(i,i)` = sum of admittances connected to bus `i`. Line
   charging susceptance (B/2) is neglected (short-line model).
2. **PQ buses** — voltage updated directly from the power-flow equation:

   ```
   V_i_new = (1/Y_ii) * (S_i / conj(V_i) - sum_j≠i Y_ij * V_j)
   ```

3. **PV buses** — reactive power `Q` is estimated each iteration from the
   current voltages, clamped to `[Qmin, Qmax]` (converting the bus to PQ if
   a limit is hit), then the voltage is updated and its magnitude is reset
   to the specified value.
4. **Acceleration** — successive over-relaxation is applied globally after
   each full sweep: `V = V_prev + alpha*(V - V_prev)`, with `alpha = 1.6`.
5. **Convergence** — stops when the largest voltage change across all
   buses drops below `tol = 1e-6`, or after `maxIter = 100` iterations.

## Sample system used for the results in this repo

3-bus system, 100 MVA base:

| Bus | Type | V (pu) | Pgen (MW) | Pload (MW) | Qload (MVAr) | Q limits (MVAr) |
|---|---|---|---|---|---|---|
| 1 | Slack | 1.05 | – | – | – | – |
| 2 | PV | 1.00 | 200 | – | – | [-100, 300] |
| 3 | PQ | – | – | 200 | 100 | – |

Lines (R, X in pu):

| From | To | R | X |
|---|---|---|---|
| 1 | 2 | 0.0200 | 0.0400 |
| 1 | 3 | 0.0100 | 0.0300 |
| 2 | 3 | 0.0125 | 0.0250 |

## Results

Converged in **26 iterations** (tolerance `1e-6`).

| Bus | Type (final) | \|V\| (pu) | Angle (deg) |
|---|---|---|---|
| 1 | Slack | 1.0500 | 0.000 |
| 2 | PQ* | 1.0000 | 0.787 |
| 3 | PQ | 1.0246 | -1.505 |

\* Bus 2 started as PV but its estimated Q went outside `[Qmin, Qmax]`
during the iterations, so the solver switched it to PQ (see "Known
caveats" below).

Run it yourself:

```bash
matlab -batch "run('gauss_seidel_loadflow_demo.m')"   # non-interactive, sample system
matlab -nosplash -r "gauss_seidel_loadflow"            # interactive, enter your own system
```

Or simply open `gauss_seidel_loadflow_demo.m` / `gauss_seidel_loadflow.m` in the MATLAB editor and click **Run**.

## Known caveats / limitations (worth mentioning in a report or viva)

1. **Units** — `Psp`/`Qsp` are divided by 100, so `Pgen`/`Qgen`/`Pload`/
   `Qload`/`Qmin`/`Qmax` must be entered in MW/MVAr on a 100 MVA base
   (not already in per-unit).
2. **PV→PQ bug** — the `pv`/`pq` index lists are built once before the
   iteration loop. When a PV bus is converted to PQ mid-iteration
   (`busdata(i,2) = 3`), those lists are never updated, so the bus keeps
   being processed by the PV branch (magnitude held fixed) for the rest of
   the run, even though `busdata` now labels it PQ. A correct fix would
   recheck `busdata(i,2)` dynamically inside the loop instead of relying on
   the precomputed lists.
3. **Final Q printout** — for the same reason, the "Reactive powers at PV
   buses" section at the end loops over the original `pv` list, so a bus
   that got converted to PQ still shows up there.
4. **No line charging** — shunt susceptance (B/2) is neglected, which is
   only valid for short lines.
5. **Relaxation placement** — over-relaxation is applied once per full
   sweep (after updating every bus), not per-bus; both approaches are used
   in practice but this is worth stating explicitly.
