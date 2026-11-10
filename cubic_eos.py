# Cubic Equation of State program (vdW, RK, SRK, PR)
# Reference: Smith, Van Ness, Abbott & Swihart (SVNA), 8th edition
# Units: T in K, P in bar, V in cm3/mol, H in J/mol, S in J/(mol K)

import math

R = 83.14       # cm3 bar / (mol K) for PVT and Poynting
R_J = 8.314     # J / (mol K) for H and S
T0 = 298.15     # reference temperature (K)
P0 = 1.0        # reference pressure (bar)

# Gas database: [Tc (K), Pc (bar), omega, A, B, C, D]
# Table B.1 (critical constants) and Table C.1 (Cp/R = A + B*T + C*T^2 + D/T^2)
gases = {
    "methane":          [190.6,  45.99,  0.012,  1.702,  9.081e-3,  -2.164e-6,  0.0],
    "ethane":           [305.3,  48.72,  0.100,  1.131, 19.225e-3,  -5.561e-6,  0.0],
    "propane":          [369.8,  42.48,  0.152,  1.213, 28.785e-3,  -8.824e-6,  0.0],
    "n-butane":         [425.1,  37.96,  0.200,  1.935, 36.915e-3, -11.402e-6,  0.0],
    "n-pentane":        [469.7,  33.70,  0.252,  2.464, 45.351e-3, -14.111e-6,  0.0],
    "n-hexane":         [507.6,  30.25,  0.301,  3.025, 53.722e-3, -16.791e-6,  0.0],
    "ethylene":         [282.3,  50.40,  0.087,  1.424, 14.394e-3,  -4.392e-6,  0.0],
    "propylene":        [365.6,  46.65,  0.140,  1.637, 22.706e-3,  -6.915e-6,  0.0],
    "benzene":          [562.2,  48.98,  0.210, -0.206, 39.064e-3, -13.301e-6,  0.0],
    "toluene":          [591.8,  41.06,  0.262,  0.290, 47.052e-3, -15.716e-6,  0.0],
    "methanol":         [512.6,  80.97,  0.564,  2.211, 12.216e-3,  -3.450e-6,  0.0],
    "ethanol":          [513.9,  61.48,  0.645,  3.518, 20.001e-3,  -6.002e-6,  0.0],
    "air":              [132.2,  37.45,  0.035,  3.355,  0.575e-3,   0.0,      -0.016e5],
    "ammonia":          [405.7, 112.80,  0.253,  3.578,  3.020e-3,   0.0,      -0.186e5],
    "carbon monoxide":  [132.9,  34.99,  0.048,  3.376,  0.557e-3,   0.0,      -0.031e5],
    "carbon dioxide":   [304.2,  73.83,  0.224,  5.457,  1.045e-3,   0.0,      -1.157e5],
    "hydrogen":         [ 33.19, 13.13, -0.216,  3.249,  0.422e-3,   0.0,       0.083e5],
    "hydrogen sulfide": [373.5,  89.63,  0.094,  3.931,  1.490e-3,   0.0,      -0.232e5],
    "nitrogen":         [126.2,  34.00,  0.038,  3.280,  0.593e-3,   0.0,       0.040e5],
    "oxygen":           [154.6,  50.43,  0.022,  3.639,  0.506e-3,   0.0,      -0.227e5],
    "sulfur dioxide":   [430.8,  78.84,  0.245,  5.699,  0.801e-3,   0.0,      -1.015e5],
    "water":            [647.1, 220.55,  0.345,  3.470,  1.450e-3,   0.0,       0.121e5],
}


def eos_params(eos):
    # Returns (u, w, Omega, Psi, sigma, eps, Zc_eos) from SVNA Table 3.1
    eos = eos.lower()
    if eos == "vdw":
        return 0.0, 0.0, 1.0 / 8.0, 27.0 / 64.0, 0.0, 0.0, 3.0 / 8.0
    if eos == "rk":
        return 1.0, 0.0, 0.08664, 0.42748, 1.0, 0.0, 1.0 / 3.0
    if eos == "srk":
        return 1.0, 0.0, 0.08664, 0.42748, 1.0, 0.0, 1.0 / 3.0
    if eos == "pr":
        s2 = math.sqrt(2.0)
        return 2.0, -1.0, 0.07780, 0.45724, 1.0 + s2, 1.0 - s2, 0.30740
    raise ValueError("Unknown EOS")


def get_kappa(eos, omega):
    eos = eos.lower()
    if eos == "srk":
        return 0.480 + 1.574 * omega - 0.176 * omega**2
    if eos == "pr":
        return 0.37464 + 1.54226 * omega - 0.26992 * omega**2
    return 0.0


def alpha_val(eos, Tr, omega):
    eos = eos.lower()
    if eos == "vdw":
        return 1.0
    if eos == "rk":
        return Tr ** (-0.5)
    k = get_kappa(eos, omega)
    return (1.0 + k * (1.0 - math.sqrt(Tr))) ** 2


def dln_alpha(eos, Tr, omega):
    # d(ln alpha) / d(ln Tr)
    eos = eos.lower()
    if eos == "vdw":
        return 0.0
    if eos == "rk":
        return -0.5
    k = get_kappa(eos, omega)
    rt = math.sqrt(Tr)
    return -k * rt / (1.0 + k * (1.0 - rt))


def calc_ab(eos, gas, T):
    Tc, Pc, omega = gases[gas][:3]
    _, _, Omega, Psi, _, _, _ = eos_params(eos)
    Tr = T / Tc
    a = Psi * alpha_val(eos, Tr, omega) * (R * Tc) ** 2 / Pc
    b = Omega * R * Tc / Pc
    return a, b


def solve_Z(eos, gas, T, P):
    if T <= 0 or P <= 0:
        raise ValueError("T and P must be positive")
    a, b = calc_ab(eos, gas, T)
    u, w, _, _, _, _, _ = eos_params(eos)

    A = a * P / (R * T) ** 2
    B = b * P / (R * T)

    al = -1.0 - B + u * B
    be = A + w * B**2 - u * B - u * B**2
    ga = -A * B - w * B**2 - w * B**3

    def f(z):
        return z**3 + al * z**2 + be * z + ga

    def df(z):
        return 3.0 * z**2 + 2.0 * al * z + be

    # Newton-Raphson starting above all roots; fallback to bisection on [B, Z_hi]
    Z_hi = max(1.0 + abs(al) + abs(be) + abs(ga), B + 1.0)
    z = Z_hi
    ok = False
    for _ in range(200):
        slope = df(z)
        if abs(slope) < 1e-14:
            break
        z_new = z - f(z) / slope
        if abs(z_new - z) < 1e-12:
            if z_new > B and abs(f(z_new)) < 1e-9:
                z = z_new
                ok = True
            break
        z = z_new

    if not ok:
        lo, hi = B, Z_hi
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if f(mid) <= 0:
                lo = mid
            else:
                hi = mid
        z = 0.5 * (lo + hi)

    z1 = z
    # Deflate (z - z1) -> quadratic z^2 + p*z + q = 0
    p = al + z1
    q = be + z1 * p
    cand = [z1]
    disc = p**2 - 4.0 * q
    if disc >= -1e-12:
        sd = math.sqrt(max(0.0, disc))
        cand.append(0.5 * (-p + sd))
        cand.append(0.5 * (-p - sd))

    # Keep physical roots Z > B (i.e. V > b) and sort ascending
    phys = sorted(r for r in cand if r > B)
    roots = []
    for r in phys:
        if not roots or abs(r - roots[-1]) > 1e-9:
            roots.append(r)
    return roots, A, B


def pressure_from_TV(eos, gas, T, V):
    if T <= 0:
        raise ValueError("T must be positive")
    a, b = calc_ab(eos, gas, T)
    if V <= b:
        raise ValueError("V must be greater than b")
    u, w, _, _, _, _, _ = eos_params(eos)
    return R * T / (V - b) - a / (V**2 + u * b * V + w * b**2)


def temperature_from_PV(eos, gas, P, V):
    if P <= 0:
        raise ValueError("P must be positive")
    Tc = gases[gas][0]
    _, b = calc_ab(eos, gas, Tc)
    if V <= b:
        raise ValueError("V must be greater than b")

    def g(t):
        return pressure_from_TV(eos, gas, t, V) - P

    T_lo, T_hi = 0.2 * Tc, 5.0 * Tc
    while g(T_lo) > 0 and T_lo > 1e-3:
        T_lo *= 0.5
    while g(T_hi) < 0 and T_hi < 1e5:
        T_hi *= 2.0
    if g(T_lo) * g(T_hi) > 0:
        raise RuntimeError("Could not bracket temperature")

    for _ in range(120):
        T_mid = 0.5 * (T_lo + T_hi)
        if T_hi - T_lo < 1e-8:
            return T_mid
        if g(T_lo) * g(T_mid) <= 0:
            T_hi = T_mid
        else:
            T_lo = T_mid
    return 0.5 * (T_lo + T_hi)


def calc_residuals(eos, gas, T, P, Z, A, B):
    Tc, _, omega = gases[gas][:3]
    _, _, _, _, sig, eps, _ = eos_params(eos)
    q = A / B
    if eos.lower() == "vdw":
        I = B / Z
    else:
        I = math.log((Z + sig * B) / (Z + eps * B)) / (sig - eps)

    d = dln_alpha(eos, T / Tc, omega)
    H_res = (Z - 1.0 + (d - 1.0) * q * I) * R_J * T
    S_res = (math.log(Z - B) + d * q * I) * R_J
    ln_phi = Z - 1.0 - math.log(Z - B) - q * I
    phi = math.exp(ln_phi)
    f = phi * P
    return H_res, S_res, ln_phi, phi, f


def icph(gas, T1, T2):
    A, B, C, D = gases[gas][3:]
    return (
        A * (T2 - T1)
        + 0.5 * B * (T2**2 - T1**2)
        + (1.0 / 3.0) * C * (T2**3 - T1**3)
        - D * (1.0 / T2 - 1.0 / T1)
    )


def icps(gas, T1, T2):
    A, B, C, D = gases[gas][3:]
    return (
        A * math.log(T2 / T1)
        + B * (T2 - T1)
        + 0.5 * C * (T2**2 - T1**2)
        - 0.5 * D * (1.0 / T2**2 - 1.0 / T1**2)
    )


def get_state(eos, gas, T, P, phase="stable"):
    roots, A, B = solve_Z(eos, gas, T, P)
    if phase == "liquid":
        Z = roots[0]
    elif phase == "vapour":
        Z = roots[-1]
    else:
        # Pick the stable root (lower ln(phi)) and discard middle root
        if len(roots) == 1:
            Z = roots[0]
        else:
            Z_L, Z_V = roots[0], roots[-1]
            lp_L = calc_residuals(eos, gas, T, P, Z_L, A, B)[2]
            lp_V = calc_residuals(eos, gas, T, P, Z_V, A, B)[2]
            Z = Z_L if lp_L <= lp_V else Z_V

    V = Z * R * T / P
    H_res, S_res, ln_phi, phi, f = calc_residuals(eos, gas, T, P, Z, A, B)
    H = R_J * icph(gas, T0, T) + H_res
    S = R_J * icps(gas, T0, T) - R_J * math.log(P / P0) + S_res
    return {
        "roots": roots, "Z": Z, "V": V, "A": A, "B": B,
        "HR": H_res, "SR": S_res, "lnphi": ln_phi, "phi": phi, "f": f,
        "H": H, "S": S,
    }


def delta_HS(eos, gas, T1, P1, T2, P2):
    s1 = get_state(eos, gas, T1, P1)
    s2 = get_state(eos, gas, T2, P2)
    dH = R_J * icph(gas, T1, T2) + s2["HR"] - s1["HR"]
    dS = R_J * icps(gas, T1, T2) - R_J * math.log(P2 / P1) + s2["SR"] - s1["SR"]
    return dH, dS


def vapour_pressure(eos, gas, T):
    Tc, Pc, omega = gases[gas][:3]
    if T <= 0 or T >= Tc:
        raise ValueError("Need 0 < T < Tc for vapour pressure")

    Zc_eos = eos_params(eos)[6]
    Vc_eos = Zc_eos * R * Tc / Pc
    P = Pc * 10.0 ** ((7.0 / 3.0) * (1.0 + omega) * (1.0 - Tc / T))
    P_lo, P_hi = None, None

    for _ in range(500):
        roots, A, B = solve_Z(eos, gas, T, P)
        if len(roots) < 2 or abs(roots[-1] - roots[0]) < 1e-8:
            V1 = roots[0] * R * T / P
            if V1 < Vc_eos:
                P_hi = P if P_hi is None else min(P_hi, P)
                P = 0.5 * (P_lo + P_hi) if P_lo is not None else 0.7 * P
            else:
                P_lo = P if P_lo is None else max(P_lo, P)
                P = 0.5 * (P_lo + P_hi) if P_hi is not None else 1.4 * P
            continue

        phi_L = calc_residuals(eos, gas, T, P, roots[0], A, B)[3]
        phi_V = calc_residuals(eos, gas, T, P, roots[-1], A, B)[3]
        P_new = P * (phi_L / phi_V)

        if P_lo is not None and P_new <= P_lo:
            P_new = 0.5 * (P + P_lo)
        if P_hi is not None and P_new >= P_hi:
            P_new = 0.5 * (P + P_hi)

        if abs(P_new - P) <= 1e-10 * P:
            P = P_new
            break
        P = P_new
    else:
        raise RuntimeError("Vapour pressure did not converge")

    liq = get_state(eos, gas, T, P, "liquid")
    vap = get_state(eos, gas, T, P, "vapour")
    return {
        "Psat": P, "VL": liq["V"], "VV": vap["V"],
        "phi": vap["phi"], "Hvap": vap["H"] - liq["H"], "Svap": vap["S"] - liq["S"],
    }


def compressed_liquid_f(eos, gas, T, P):
    sat = vapour_pressure(eos, gas, T)
    poynting = math.exp(sat["VL"] * (P - sat["Psat"]) / (R * T))
    f_poy = sat["Psat"] * sat["phi"] * poynting
    f_eos = get_state(eos, gas, T, P, "liquid")["f"]
    return {
        "f_poynting": f_poy, "poynting": poynting,
        "Psat": sat["Psat"], "phi_sat": sat["phi"], "VL": sat["VL"], "f_eos": f_eos,
    }


def check_density_integrals(eos, gas, T, P, Z, N=20000):
    V = Z * R * T / P
    rho_max = 1.0 / V
    dr = rho_max / N
    u, w, _, _, _, _, _ = eos_params(eos)
    dT = max(1e-4, 1e-5 * T)

    def Z_rho(t, rho):
        a_t, b_t = calc_ab(eos, gas, t)
        return 1.0 / (1.0 - b_t * rho) - (a_t * rho) / (
            R * t * (1.0 + u * b_t * rho + w * (b_t * rho) ** 2)
        )

    def eval_pt(rho):
        if rho == 0.0:
            a0, b0 = calc_ab(eos, gas, T)
            ap, _ = calc_ab(eos, gas, T + dT)
            am, _ = calc_ab(eos, gas, T - dT)
            t1 = b0 - a0 / (R * T)
            t2 = -(ap / (R * (T + dT)) - am / (R * (T - dT))) / (2.0 * dT)
            return t1, t2
        z0 = Z_rho(T, rho)
        zp = Z_rho(T + dT, rho)
        zm = Z_rho(T - dT, rho)
        return (z0 - 1.0) / rho, (zp - zm) / (2.0 * dT * rho)

    a1, a2 = eval_pt(0.0)
    b1, b2 = eval_pt(rho_max)
    s1 = 0.5 * (a1 + b1)
    s2 = 0.5 * (a2 + b2)
    for i in range(1, N):
        v1, v2 = eval_pt(i * dr)
        s1 += v1
        s2 += v2

    I1 = s1 * dr
    I2 = s2 * dr
    HR_num = (Z - 1.0 - T * I2) * R_J * T
    SR_num = (math.log(Z) - T * I2 - I1) * R_J
    lp_num = Z - 1.0 - math.log(Z) + I1
    return HR_num, SR_num, lp_num


def run_validations():
    print("\n1) SVNA Example 3.9: n-butane, RK, 350 K, 9.4573 bar")
    T, P = 350.0, 9.4573
    roots, A, B = solve_Z("rk", "n-butane", T, P)
    VL, VV = roots[0] * R * T / P, roots[-1] * R * T / P
    print(f"   A = {A:.6f}, B = {B:.6f}, roots = {[round(r, 6) for r in roots]}")
    print(f"   VL = {VL:.2f} cm3/mol (Ex 3.9 ~ 133.3), VV = {VV:.2f} cm3/mol (Ex 3.9 ~ 2555)")
    print(f"   Round-trip P(T,VV) = {pressure_from_TV('rk', 'n-butane', T, VV):.6f} bar, "
          f"T(P,VV) = {temperature_from_PV('rk', 'n-butane', P, VV):.4f} K")

    print("\n2) Ideal-gas limit Z -> 1 as P -> 0 (n-butane, 350 K)")
    print(f"{'P (bar)':>10} | {'vdW':>10} | {'RK':>10} | {'SRK':>10} | {'PR':>10}")
    print("-" * 62)
    for Pt in (10.0, 1.0, 0.1, 0.01, 1e-4):
        row = [f"{solve_Z(e, 'n-butane', 350.0, Pt)[0][-1]:10.6f}" for e in ("vdw", "rk", "srk", "pr")]
        print(f"{Pt:10.4g} | " + " | ".join(row))

    print("\n3) Consistency (HR/RT - SR/R == ln(phi)) and Trapezoidal Integral Check (N=20,000)")
    print(f"{'State':<28} | {'EOS':^4} | {'|Diff|':>8} | {'HR (cl / num)':^19} | {'SR (cl / num)':^17} | {'ln(phi) (cl / num)':^19}")
    print("-" * 109)
    cases = [
        ("n-butane", "vdw", 500.0, 50.0, "stable"),
        ("n-butane", "rk",  500.0, 50.0, "stable"),
        ("n-butane", "srk", 500.0, 50.0, "stable"),
        ("n-butane", "pr",  500.0, 50.0, "stable"),
        ("carbon dioxide", "srk", 320.0, 100.0, "stable"),
        ("methane", "rk", 200.0, 50.0, "stable"),
        ("n-butane", "rk", 350.0, 9.4573, "liquid"),
    ]
    for gas, eos, Tc_v, Pc_v, ph in cases:
        st = get_state(eos, gas, Tc_v, Pc_v, ph)
        diff = abs(st["HR"] / (R_J * Tc_v) - st["SR"] / R_J - st["lnphi"])
        hn, sn, ln = check_density_integrals(eos, gas, Tc_v, Pc_v, st["Z"], 20000)
        lbl = f"{gas} ({Tc_v:g}K,{Pc_v:g}b,{ph[0]})"
        print(f"{lbl:<28} | {eos.upper():^4} | {diff:8.2e} | {st['HR']:8.2f}/{hn:8.2f} | "
              f"{st['SR']:7.4f}/{sn:7.4f} | {st['lnphi']:8.5f}/{ln:8.5f}")

    print("\n4) Vapour pressure & Poynting test (n-butane, PR, 350 K)")
    sat = vapour_pressure("pr", "n-butane", 350.0)
    liq = compressed_liquid_f("pr", "n-butane", 350.0, 50.0)
    print(f"   Psat = {sat['Psat']:.4f} bar, VL = {sat['VL']:.2f}, VV = {sat['VV']:.2f} cm3/mol, phi = {sat['phi']:.5f}")
    print(f"   At 50 bar: f_Poynting = {liq['f_poynting']:.4f} bar, f_EOS = {liq['f_eos']:.4f} bar")


def plot_Z_vs_P(gas="n-butane", T=500.0, P_max=150.0, filename="fig_Z_vs_P.png"):
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib is not installed.")
        return

    ps = [0.5 + (P_max - 0.5) * i / 149.0 for i in range(150)]
    styles = [("vdw", "vdW", "--"), ("rk", "RK", "-."), ("srk", "SRK", ":"), ("pr", "PR", "-")]

    plt.figure(figsize=(6.8, 4.2), dpi=200)
    for eos_name, label, ls in styles:
        zs = [get_state(eos_name, gas, T, p)["Z"] for p in ps]
        plt.plot(ps, zs, linestyle=ls, linewidth=2.0, label=label)

    plt.axhline(1.0, color="black", linestyle=":", linewidth=1.0, alpha=0.6, label="Ideal gas (Z = 1)")
    plt.xlabel("Pressure, P (bar)")
    plt.ylabel("Compressibility Factor, Z")
    plt.title(f"Compressibility Factor Z vs P for {gas} at T = {T:g} K")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()
    print(f"Saved plot to {filename}")


def pick_gas_eos():
    gas = input("\nEnter gas name: ").strip().lower()
    while gas not in gases:
        gas = input("Not in list, try again: ").strip().lower()
    eos = input("Enter EOS (vdw / rk / srk / pr): ").strip().lower()
    while eos not in ("vdw", "rk", "srk", "pr"):
        eos = input("Enter vdw, rk, srk, or pr: ").strip().lower()
    return gas, eos


def main():
    print("Cubic EOS Calculator (vdW, RK, SRK, PR)")
    print("Available gases:", ", ".join(gases.keys()))
    gas, eos = pick_gas_eos()

    while True:
        print(f"\n[{gas}, {eos.upper()}]")
        print("  1  PVT, H, S and fugacity at (T, P)")
        print("  2  Pressure from (T, V)")
        print("  3  Temperature from (P, V)")
        print("  4  Vapour pressure at T")
        print("  5  Compressed-liquid fugacity (Poynting)")
        print("  6  Change in H and S between two states")
        print("  7  Run validation checks")
        print("  8  Plot Z vs P for all 4 EOS")
        print("  9  Change gas / EOS")
        print("  0  Exit")
        ch = input("Choice: ").strip()

        if ch == "0":
            break
        elif ch == "1":
            try:
                T = float(input("T (K): "))
                P = float(input("P (bar): "))
                if T < 298.15:
                    print("Note: T < 298.15 K extrapolates Table C.1 Cp/R fit.")
                s = get_state(eos, gas, T, P)
                print(f"A = {s['A']:.6f}, B = {s['B']:.6f}")
                print("Roots (Z > B):", [round(r, 6) for r in s["roots"]])
                if len(s["roots"]) == 3:
                    print(f"Middle root Z_M = {s['roots'][1]:.6f} (discarded)")
                print(f"Selected Z = {s['Z']:.6f}, V = {s['V']:.3f} cm3/mol")
                print(f"HR = {s['HR']:.2f} J/mol, SR = {s['SR']:.4f} J/(mol K)")
                print(f"ln(phi) = {s['lnphi']:.6f}, phi = {s['phi']:.6f}, f = {s['f']:.4f} bar")
                print(f"H = {s['H']:.2f} J/mol, S = {s['S']:.4f} J/(mol K)  [ref: ig 298.15 K, 1 bar]")
            except ValueError as e:
                print("Error:", e)
        elif ch == "2":
            try:
                T = float(input("T (K): "))
                V = float(input("V (cm3/mol): "))
                print(f"P = {pressure_from_TV(eos, gas, T, V):.5f} bar")
            except ValueError as e:
                print("Error:", e)
        elif ch == "3":
            try:
                P = float(input("P (bar): "))
                V = float(input("V (cm3/mol): "))
                print(f"T = {temperature_from_PV(eos, gas, P, V):.4f} K")
            except Exception as e:
                print("Error:", e)
        elif ch == "4":
            try:
                T = float(input("T (K): "))
                sat = vapour_pressure(eos, gas, T)
                print(f"Psat = {sat['Psat']:.4f} bar, phi = {sat['phi']:.5f}")
                print(f"VL = {sat['VL']:.3f} cm3/mol, VV = {sat['VV']:.2f} cm3/mol")
                print(f"dHvap = {sat['Hvap']:.1f} J/mol, dSvap = {sat['Svap']:.4f} J/(mol K)")
            except Exception as e:
                print("Error:", e)
        elif ch == "5":
            try:
                T = float(input("T (K): "))
                P = float(input("P (bar): "))
                res = compressed_liquid_f(eos, gas, T, P)
                if P < res["Psat"]:
                    print(f"Note: P < Psat ({res['Psat']:.4f} bar), stable phase is vapour.")
                print(f"Psat = {res['Psat']:.4f} bar, Poynting = {res['poynting']:.5f}")
                print(f"f (Poynting) = {res['f_poynting']:.4f} bar, f (EOS liquid) = {res['f_eos']:.4f} bar")
            except Exception as e:
                print("Error:", e)
        elif ch == "6":
            try:
                T1 = float(input("T1 (K): "))
                P1 = float(input("P1 (bar): "))
                T2 = float(input("T2 (K): "))
                P2 = float(input("P2 (bar): "))
                dH, dS = delta_HS(eos, gas, T1, P1, T2, P2)
                print(f"Delta H = {dH:.2f} J/mol, Delta S = {dS:.4f} J/(mol K)")
            except ValueError as e:
                print("Error:", e)
        elif ch == "7":
            run_validations()
        elif ch == "8":
            try:
                T = float(input("T (K): "))
                P_max = float(input("P_max (bar): "))
                plot_Z_vs_P(gas, T, P_max)
            except ValueError as e:
                print("Error:", e)
        elif ch == "9":
            gas, eos = pick_gas_eos()
        else:
            print("Invalid choice.")


if __name__ == "__main__":
    main()
