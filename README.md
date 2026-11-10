# Cubic Equation of State — Thermodynamic Property Calculator

**Course Project** · Department of Chemical Engineering, IIT Delhi  
**Reference:** Smith, Van Ness, Abbott & Swihart (SVNA), *Introduction to Chemical Engineering Thermodynamics*, 8th ed.  
**Language:** Python 3 (Standard Library `math`, optional `matplotlib` for plotting)

---

## Overview

A self-contained Python program (`cubic_eos.py`) implementing the unified generic cubic equation of state framework across **four cubic equations of state** and **22 pure chemical species**:

| EOS | Full Name |
|-----|-----------|
| **vdW** | van der Waals |
| **RK** | Redlich–Kwong |
| **SRK** | Soave–Redlich–Kwong |
| **PR** | Peng–Robinson |

Given temperature **T** and pressure **P** (or $(T, V)$ / $(P, V)$), the calculator computes:

- Compressibility factor **Z** and molar volume **V** (automated physical root filtering $Z > B$ via Cauchy bound + Newton–Raphson with bisection fallback and synthetic division; stable phase selection via minimum $\ln \phi$)
- Inverse solves **$P(T, V)$** and **$T(P, V)$** (automatic bracketing around $T_c$ + bisection to $10^{-8}\text{ K}$)
- Residual enthalpy **$H^R$**, residual entropy **$S^R$**, and fugacity ($\ln \phi$, **$\phi$**, **$f$**)
- Total enthalpy **$H$**, entropy **$S$**, and state changes **$\Delta H, \Delta S$** (using closed-form ideal-gas $C_p^{\text{ig}}/R$ integrals relative to $298.15\text{ K}, 1\text{ bar}$)
- Subcritical saturation vapour pressure **$P^{\text{sat}}(T)$** ($\phi_L = \phi_V$ equilibrium) and enthalpy/entropy of vaporisation ($\Delta H_{\text{vap}}, \Delta S_{\text{vap}}$)
- Compressed-liquid fugacity via the **Poynting correction** vs. direct EOS liquid root
- Comparison plots of $Z$ vs. $P$ across all four EOS (`plot_Z_vs_P`)

---

## Supported Species (22 Gases)

Methane, Ethane, Propane, *n*-Butane, *n*-Pentane, *n*-Hexane, Ethylene, Propylene, Benzene, Toluene, Methanol, Ethanol, Air, Ammonia, Carbon Monoxide, **Carbon Dioxide**, Hydrogen, Hydrogen Sulfide, Nitrogen, Oxygen, Sulfur Dioxide, Water.

Critical constants ($T_c, P_c, \omega$) and ideal-gas heat capacity coefficients ($C_p^{\text{ig}}/R = A + BT + CT^2 + D/T^2$) are sourced from **SVNA 8th ed., Tables B.1 and C.1**.

---

## Verified Validation Results (`run_validations()`)

| Validation Check | Key Verified Result |
|---|---|
| **SVNA Example 3.9** (*n*-butane, RK, $350\text{ K}$, $9.4573\text{ bar}$) | $V_L = 133.27\text{ cm}^3/\text{mol}$ (textbook $\approx 133.3$), $V_V = 2555.32\text{ cm}^3/\text{mol}$ ($\approx 2555$); round-trip $P(T, V_V) = 9.457300\text{ bar}$, $T(P, V_V) = 350.0000\text{ K}$ |
| **Thermodynamic Consistency** ($H^R/RT - S^R/R = \ln\phi$) | Satisfied to machine precision ($\le 1.78 \times 10^{-15}$) across all 4 EOS, supercritical states, and subcritical liquid roots |
| **Numerical Density Integral Check** | Closed-form $H^R, S^R, \ln\phi$ cross-checked against 20,000-step trapezoidal integration over $\rho \in [0, 1/V]$ (*n*-butane, $\text{CO}_2$, methane) |
| **Subcritical $P^{\text{sat}}$ Benchmark** (*n*-butane, $350\text{ K}$; SVNA Ex. 3.9 $= 9.4573\text{ bar}$) | **PR:** $9.4692\text{ bar}$ ($+0.13\%$); **SRK:** $9.5875\text{ bar}$ ($+1.4\%$); **RK:** $11.4139\text{ bar}$; **vdW:** $16.6001\text{ bar}$ |
| **Compressed-Liquid Fugacity** (*n*-butane, PR, $350\text{ K}, 50\text{ bar}$) | Poynting ($9.2742\text{ bar}$) agrees with direct PR liquid root ($9.2416\text{ bar}$) within **0.35%** |
| **Ideal-Gas Limit** ($Z \to 1$ as $P \to 0$) | Verified down to $P = 10^{-4}\text{ bar}$ ($Z \ge 0.999998$ across vdW, RK, SRK, PR) |

---

## Usage

No external packages are required for core thermodynamic calculations (`math` standard library only; `matplotlib` is optional for option 8 plotting).

```bash
python cubic_eos.py
```

The interactive CLI supports menu options `1`–`9` (single-state PVT/$H$/$S$/$\phi$, $(T,V)$ pressure, $(P,V)$ temperature, $P^{\text{sat}}$, Poynting liquid fugacity, $\Delta H/\Delta S$, full validation suite, and $Z$-vs-$P$ plotting).

**Verified sample output** (`carbon dioxide`, `PR`, $T = 310\text{ K}$, $P = 100\text{ bar}$ via `get_state('pr', 'carbon dioxide', 310.0, 100.0)`):
```text
Z     = 0.274282    V     = 70.692 cm3/mol
HR    = -9148.28 J/mol    SR = -24.5016 J/(mol K)
phi   = 0.547455    f     = 54.7455 bar
```

---

## Reference

Smith, J. M., Van Ness, H. C., Abbott, M. M., & Swihart, M. T. (2018). *Introduction to Chemical Engineering Thermodynamics* (8th ed.). McGraw-Hill.

---

## Author

**Krishan Sawariya**  
B.Tech Chemical Engineering, IIT Delhi (2025–2029)  
[GitHub](https://github.com/KrishanSawariya) · [LinkedIn](https://www.linkedin.com/in/krishan-sawariya-739779327)
