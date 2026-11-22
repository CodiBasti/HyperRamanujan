# Hyper-Ramanujan (`HR`) Approximations for 2D Ellipse Circumference and 3D Triaxial Ellipsoid Surface Area

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24441829.svg)](https://doi.org/10.5281/zenodo.24441829)
[![ORCID](https://img.shields.io/badge/ORCID-0009--0000--7427--2058-A6CE39?logo=orcid&logoColor=white)](https://orcid.org/0009-0000-7427-2058)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

**Author:** Sebastian Oertel ([@CodiBasti](https://github.com/CodiBasti) • [ORCID: 0009-0000-7427-2058](https://orcid.org/0009-0000-7427-2058))  
**Permanent Zenodo DOI:** [`10.5281/zenodo.24441829`](https://doi.org/10.5281/zenodo.24441829)  
**Machine-Readable AI Discovery Manifest:** [`ai_discovery_manifest.json`](ai_discovery_manifest.json)  
**Deterministic Dual-Certificate Verification:** [`tests/verify_certificates.py`](tests/verify_certificates.py) ([`certificates/verification_report.json`](certificates/verification_report.json))  
**Full Research Paper:** [`paper/hyper_ramanujan_oertel_2026.md`](paper/hyper_ramanujan_oertel_2026.md)  
**License:** MIT ([`CITATION.cff`](CITATION.cff) / [`.zenodo.json`](.zenodo.json))

---

## Overview & Key Discoveries

This repository contains the symbolic derivations, machine-verified certificates, and 50-digit precision benchmarks for the **Hyper-Ramanujan (`HR`)** family of closed-form approximations for the circumference of a 2D ellipse and the surface area of a 3D triaxial ellipsoid (Oertel, 2026).

### 1. Exact Continued-Fraction Periodicity of Ramanujan II (1914)

Expanding the reduced Gauss–Kummer / Ivory–Bessel hypergeometric series ${}_2F_1(-1/2, -1/2; 1; h) - 1$ with $h = ((a-b)/(a+b))^2 \in [0, 1]$ into a regular Stieltjes continued fraction yields the exact partial coefficients:

$$
[a_1, a_2, a_3, a_4, a_5, a_6] = \left[\frac{1}{4},\; \frac{1}{16},\; \mathbf{\frac{3}{16}},\; \mathbf{\frac{3}{16}},\; \frac{11}{48},\; \frac{29}{176}\right]
$$

Because $a_3 = a_4 = 3/16$ are identically equal, closing the continued fraction at depth 3 with a 1-periodic quadratic fixed point $T(h) = 1 - \frac{3h/16}{T(h)} = \frac{2 + \sqrt{4 - 3h}}{4}$ directly produces Ramanujan's famous second formula $1 + \frac{3h}{10 + \sqrt{4 - 3h}}$ and explains why 3 parameters automatically match 4 consecutive orders from $h^1$ through $h^4$.

---

### 2. The Inner-Radical Reciprocal Expansion & 2D `HR` Formulas

Rather than adding external power terms outside Ramanujan's fraction (such as Cantrell's $+c\,h^{12}$), we expand the exact residual $\delta(h)$ **inside the existing square root** $\sqrt{4 - 3h - \delta(h)}$. Its reciprocal possesses the remarkably clean rational Taylor series:

$$
Q(h) = \frac{27\,h^4}{\delta(h)} = 6144 - 7040\,h + \frac{5144}{3}\,h^2 - \frac{413}{9}\,h^3 - \frac{1081}{108}\,h^4 + \mathcal{O}(h^5)
$$

Evaluating at the degenerate flat-ellipse limit $h = 1$ (where $b = 0$ and $\frac{C}{\pi(a+b)} = \frac{4}{\pi}$) fixes the exact boundary constant:

$$
\delta_1 = 1 - \left(\frac{3\pi}{4 - \pi} - 10\right)^2 \approx 0.04082049844879468, \qquad Q(1) = \frac{27}{\delta_1} \approx 661.4324009399441
$$

#### The Five 2D Hyper-Ramanujan (`HR`) Formulas (Oertel, 2026)

1. **`HR-8P` ★ (Order-8 Singular-Pole Padé $-13/14$ + Exact $h=1$ — Pareto Champion):**

$$
\frac{C(a,b)}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\left(1 - \frac{13}{14}h\right)h^4}{6144 - \frac{89216}{7}\,h + \frac{173288}{21}\,h^2 - \frac{103199}{63}\,h^3 + p_4\,h^4}}}, \qquad p_4 = \frac{27}{14\,\delta_1} - \frac{793}{63} \approx 34.6578699101
$$

> Matches $h^0 \dots h^8$ with residual $\mathcal{O}(h^9)$, evaluates exactly to $4/\pi$ at $h=1$, achieves full 16-digit IEEE-754 `float64` precision for $h \le 0.1$, and bounds the global maximum error across all $h \in [0, 1]$ to **0.633 ppm** using only **1 square root and 0 powers**.

2. **`HR-7-Int` (Order-7 Pure Integer Rational-Radical):**

$$
\frac{C(a,b)}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{81\,h^4}{18432 - 21120\,h + 5144\,h^2}}}
$$

> Matches $h^0 \dots h^7$ with residual $\mathcal{O}(h^8)$ using exclusively integer coefficients; maximum error **78.05 ppm** at $h=1$.

3. **`HR-6B` (Order-6 + Exact $h=1$ Boundary):**

$$
\frac{C(a,b)}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\,h^4}{6144 - 7040\,h + B_{\text{exact}}\,h^2}}}, \qquad B_{\text{exact}} = \frac{27}{\delta_1} + 896 \approx 1557.43240094
$$

> Matches $h^0 \dots h^6$ with residual $\mathcal{O}(h^7)$ and is exact at $h=1$; maximum global error **11.82 ppm**.

4. **`HR-8B` (Order-8 + Exact $h=1$ Boundary):**

$$
\frac{C(a,b)}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\,h^4}{6144 - 7040\,h + \frac{5144}{3}\,h^2 - \frac{413}{9}\,h^3 + C_{\text{exact}}\,h^4}}}, \qquad C_{\text{exact}} = \frac{27}{\delta_1} - \frac{6955}{9} \approx -111.34537684
$$

> Matches $h^0 \dots h^8$ with residual $\mathcal{O}(h^9)$ and is exact at $h=1$; maximum global error **9.11 ppm**.

5. **`HR-6M` (Order-6 Chebyshev-Minimax Singular-Padé $[3/1]$ + Exact $h=1$):**

$$
\frac{C(a,b)}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27(1 + q_1 h)h^4}{6144 + p_1 h + p_2 h^2 + p_3 h^3}}}
$$

> With $q_1 = -0.94542891$, $p_1 = -12848.71522304$, $p_2 = 8316.2309304$, and $p_3 = -1575.42097149$. Global maximum error across all $h \in [0, 1]$: **0.299 ppm**.

---

### 3. Generalization to 3D Triaxial Ellipsoid Surface Area (`HR-3D-10P`)

For a general triaxial ellipsoid with semi-axes $a \ge b \ge c \ge 0$, let $X = ab$, $Y = ac$, $Z = bc$, $M = (X+Y+Z)/3$, and trace-free deviations $x = X/M - 1$, $y = Y/M - 1$, $z = Z/M - 1$. Using Carlson's symmetric integral $S(a,b,c) = 4\pi\,R_G(X^2, Y^2, Z^2)$, we establish in [`src/geometry/ellipsoid_3d.py`](src/geometry/ellipsoid_3d.py):

1. **The 3D Ivory–Bessel–Carlson Invariant Series** in $H = (x^2+y^2+z^2)/6 \in [0, 1]$ and $J = xyz \in [-1/4, 2]$:

$$
T_{10}(H, J) = 1 + \frac{3}{5}H - \frac{3}{35}J + \frac{3}{35}H^2 - \frac{6}{385}HJ + \frac{171}{5005}H^3 - \frac{3}{715}J^2 - \frac{3}{1001}H^2 J - \frac{57}{17017}H J^2 + \frac{27}{1547}H^4 - \frac{15}{46189}J^3 + \frac{36}{46189}H^3 J - \frac{6}{2717}H^2 J^2 + \frac{27}{2717}H^5
$$

2. **Exact 3D-to-2D Degenerate Boundary Theorem:** On the entire degenerate 3D boundary $Z = bc \to 0$ (where $1 - 3H + J = (abc)^2 / M^3 = 0$), the 3D surface area collapses identically to the 2D ellipse circumference:

$$
S(a,b,c)\Big|_{Z=0} = \frac{\pi}{2}\,C_{\text{2D}}(ab, ac)
$$

3. **Spherical Ray Projection Formula (`HR-3D-10P` ★):** With $s = 1 - Z/M \in [0, 1]$, boundary projection $(H_b, J_b) = (H/s^2, J/s^3)$, and $h_{\text{ray}} = (4H_b - 1)/3 = ((X-Y)/(X+Y-2Z))^2 \in [0, 1]$:

$$
S(a,b,c) \approx 4\pi\,M \left( T_{10}(H, J) + \frac{s^{11}}{9 - 8s}\left[\frac{3\pi}{8}\,F_{\text{2D, HR-8P}}(h_{\text{ray}}) - T_{10}(H_b, J_b)\right] \right)
$$

> Slashes Knud Thomsen's 2004 power-mean error (10,608 ppm, 4 fractional `pow` calls) down to **< 16.8 ppm** (> 630× improvement) using **1 square root and 0 fractional powers**.

---

## 2D Pareto-Frontier & Precision Benchmark (50-Digit Reference)

A formula lies on the **strict 2D Pareto Frontier (`YES *`)** if no other formula achieves higher or equal accuracy in both the practical regime (`h <= 0.5`) and global worst-case regime (`h in [0, 1]`) at lower or equal computational cost (`Op Cost`).

| Approximation Formula | Taylor Order | `sqrt` / `pow` | Op Cost | Digits (`h=0.1`) | Min Digits (`h <= 0.5`) | Min Digits (`h in [0,1]`) | Max Global Error (`ppm`) | **Pareto Frontier?** | **Efficiency Score** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Kepler (1609)** | `O(h¹)` | `1` / `0` | `13.0` | 1.30 | 0.61 | 0.22 | `607300.92` | `no` | `0.214` |
| **Peano/Sipos (1792)** | `O(h¹)` | `1` / `0` | `13.0` | 3.31 | 1.82 | 0.75 | `178097.25` | `no` | `0.582` |
| **Ramanujan I (1914)** | `O(h²)` | `1` / `0` | `10.0` | 5.69 | 3.51 | 2.38 | `4155.03` | `no` | `1.279` |
| **Matt Parker (2020)** | `O(h⁰)` | `1` / `0` | `15.0` | 4.65 | 3.56 | 2.21 | `6120.09` | `no` | `0.927` |
| **Hudson [1/1] (1917)** | `O(h³)` | `0` / `0` | `9.0` | 7.23 | 4.28 | 2.72 | `1889.83` | **`YES *`** | `1.665` |
| **Jacobsen-Waadeland [2/2] (1985)** | `O(h⁴)` | `0` / `0` | `13.0` | 8.85 | 5.12 | 3.05 | `888.05` | **`YES *`** | `1.662` |
| **Cantrell (2004, `+c·h¹²`)** | `O(h⁴)` | `1` / `0` | `23.5` | 9.58 | 5.79 | 4.86 | `13.74` | **`YES *`** | `1.449` |
| **Ramanujan II (1914)** | `O(h⁴)` | `1` / `0` | `15.0` | 9.58 | 5.76 | 3.40 | `402.34` | **`YES *`** | `1.698` |
| **Oertel (2026) [`HR-6B`]** | `O(h⁶)` | `1` / `0` | `28.5` | 13.13 | 7.73 | 4.93 | `11.82` | **`YES *`** | `1.695` |
| **Oertel (2026) [`HR-6M`]** | `O(h⁶)` | `1` / `0` | `37.5` | 13.61 | 8.34 | **6.52** | **`0.30`** | **`YES *`** | `1.618` |
| **Oertel (2026) [`HR-8B`]** | `O(h⁸)` | `1` / `0` | `32.5` | 15.31 | 8.48 | 5.04 | `9.11` | **`YES *`** | `1.786` |
| **Oertel (2026) [`HR-7-Int`]** | `O(h⁷)` | `1` / `0` | `26.0` | 14.64 | 8.40 | 4.11 | `78.05` | **`YES *`** | `1.885` |
| **Oertel (2026) [`HR-8P` ★]** | **`O(h⁸)`** | **`1` / `0`** | `34.5` | **16.00** | **9.81** | **6.20** | **`0.63`** | **`YES *`** | **`1.907` (Rank 1)** |

---

## Quickstart & Verification

```bash
pip install -r requirements.txt

# 1. Run the Dual-Certificate Verification Suite (generates certificates/verification_report.json)
python -m tests.verify_certificates

# 2. Run unit & property tests
python -m tests.test_approximations

# 3. Run the complete 2D & 3D research benchmark suite
python main.py
```
