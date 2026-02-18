# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Running the Program

```bash
pip install -r requirements.txt
python evaporator.py
```

There are no tests or linting tools configured. The Excel export block at the bottom of `evaporator.py` is commented out — uncomment it to write results to `evaporator_results.xlsx`.

## Project Overview

A single-script simulation tool for **sugar industry multiple-effect evaporators**, implementing Harold Birkett's method (single-set solver) extended by Cody Eaves to handle multiple parallel sets. It is **not** a general-purpose evaporator tool — the thermodynamic correlations are calibrated specifically for sugar juice.

All units are US customary: tons/hr (tph), °F, psia, ft², BTU/lb.

## Architecture — `evaporator.py`

The file is structured in three logical sections:

### 1. Input Configuration (lines 22–73)
All process inputs are hardcoded module-level variables:
- Feed juice conditions, desired syrup brix, exhaust steam pressures, last-effect vacuum
- Per-set lists: heating surface (`hs_list`), vapor bleeds (`vbleed_list`)
- Boolean flags (`set_1_online`, etc.) control which sets and the pre-evaporator are active

### 2. Thermodynamic Utility Functions (lines 83–257)
Pure functions used throughout the solver:
- `sat_steam_temp(p_psia)` / `get_latent_heat(p_psia)` — fast polynomial approximations valid only for **1–60 psia**
- `bpe_brix()`, `bpe_head()`, `bpe_total()` — boiling point elevation for sugar juice
- `get_cp(brix)` — specific heat
- `shortcut_evap()` — mass-balance shortcut to seed iterative solvers
- `u_dessin()` / `u_calc()` — Dessin and direct heat-transfer coefficient methods

### 3. Solvers and Execution (lines 264–776)

**`solve_set()`** — iterative solver for a single evaporator set (3–5 effects):
- Outer loop (50 iterations): adjusts intermediate pressure profile using the ratio of U_calc to U_dessin per effect, nudging pressures toward a uniform U-ratio
- Middle loop (5 trials per outer iteration): recalculates brix, BPE, temperature, and cp profiles
- Inner `while` loop: converges exhaust steam flow (`exh_in`) to match required total evaporation

**`pre_evap()`** — simplified 20-iteration solver for the single-effect pre-evaporator; vapor bleed is fixed (not solved for), and the loop converges vapor pressure

**Multi-set juice distribution optimizer** (lines 668–745):
- `objective_function()` calls `solve_set()` for each active set at trial juice fractions and returns the pairwise differences in mean U-ratio
- `scipy.optimize.fsolve` drives these differences to zero, finding optimal juice splits
- Results are stored in `df_set_1_opt`, `df_set_2_opt`, `df_set_3_opt` (pandas DataFrames) and `df_pre_3` (pandas Series)

### Key Constraints to Keep in Mind
- Steam property polynomials are only accurate for 1–60 psia; inputs outside this range silently produce wrong results
- `v_bld_list` must always have exactly 3 elements (for effects 1–3); the code appends zeros for 4th/5th effects internally
- The pre-evaporator result feeds brix and juice temperature into all downstream sets when `pre_3_online = True`
