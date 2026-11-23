# Closed-Form Radical Approximations for the Ellipse Circumference and Triaxial Ellipsoid Surface Area: Continued Fraction Periodicity of Ramanujan II and Singular-Pole Hyper-Ramanujan (`HR`) Extensions

**Author:** Sebastian Oertel ([@CodiBasti](https://github.com/CodiBasti) • [ORCID: 0009-0000-7427-2058](https://orcid.org/0009-0000-7427-2058))  
**Preprint DOI:** [`10.5281/zenodo.24462362`](https://doi.org/10.5281/zenodo.24462362) • **Code DOI:** [`10.5281/zenodo.24441829`](https://doi.org/10.5281/zenodo.24441829)  
**Repository:** `https://github.com/CodiBasti/HyperRamanujan`  
**Date:** November 2026  

---

## Abstract

We investigate the analytical origin of Srinivasa Ramanujan's famous 1914 approximation for the circumference of an ellipse,

$$
C(a,b) \approx \pi(a+b)\left(1 + \frac{3h}{10 + \sqrt{4 - 3h}}\right), \qquad h = \left(\frac{a-b}{a+b}\right)^2
$$

and resolve the long-standing question of why a three-parameter radical ansatz matches four consecutive orders ($h^1$ through $h^4$) of the Ivory–Bessel hypergeometric expansion ${}_2F_1(-1/2, -1/2; 1; h)$. By expanding ${}_2F_1(-1/2, -1/2; 1; h) - 1$ into a regular Stieltjes continued fraction, we prove that the third and fourth partial numerators are **identically equal** ($a_3 = a_4 = 3/16$). Closing the continued fraction tail at level 3 with a 1-periodic quadratic fixed point $T(h) = 1 - \frac{3h/16}{T(h)}$ yields Ramanujan's formula directly and explains the automatic cancellation of the $h^4$ error defect.

Building on this structure, we analyze the exact residual $\delta(h)$ **inside the square root** $\sqrt{4 - 3h - \delta(h)}$ and show that its reciprocal possesses the remarkably clean rational expansion

$$
\frac{27\,h^4}{\delta(h)} = 6144 - 7040\,h + \frac{5144}{3}\,h^2 - \frac{413}{9}\,h^3 - \frac{1081}{108}\,h^4 + \mathcal{O}(h^5)
$$

From this expansion we derive five novel single-radical Hyper-Ramanujan approximations (**`HR-7-Int`**, **`HR-6B`**, **`HR-8B`**, **`HR-8P`**, and **`HR-6M`**), culminating in the singular-pole Padé radical **`HR-8P`** with pole factor $1 - \frac{13}{14}h$ that matches the Taylor series up to $\mathcal{O}(h^9)$, evaluates **identically to the exact degenerate limit $4/\pi$ at $h=1$**, reaches full IEEE-754 `float64` machine precision ($1.15 \times 10^{-17}$) at moderate eccentricities, and slashes the global maximum error across all $h \in [0, 1]$ from 402.34 ppm (Ramanujan II) down to **0.63 ppm** (**`HR-8P`**) and **0.30 ppm** (**`HR-6M`**) using only a single square root. Finally, we generalize the method to the 3D surface area of a triaxial ellipsoid (**`HR-3D-10P`**), proving that the entire 3D degenerate boundary reduces identically to $\frac{\pi}{2}C_{\text{2D}}(ab, ac)$ and slashing Knud Thomsen's 2004 power-mean error (10,608 ppm) down to **< 16.8 ppm** without fractional powers.

---

## 1. Exact Formulation and the Degenerate Limit $4/\pi$ at $h = 1$

Let an ellipse have semi-major axis $a > 0$ and semi-minor axis $b \in [0, a]$, with symmetric shape parameter

$$
h = \left(\frac{a-b}{a+b}\right)^2 \in [0, 1]
$$

Using Gauss's quadratic transformation of the complete elliptic integral of the second kind $E(e^2)$, James Ivory (1796) and Friedrich Wilhelm Bessel (1825) established the rapidly converging series representation

$$
\frac{C(a,b)}{\pi(a+b)} = {}_2F_1\!\left(-\frac{1}{2}, -\frac{1}{2}; 1; h\right) = \sum_{n=0}^{\infty} \binom{1/2}{n}^2 h^n = 1 + \frac{1}{4}h + \frac{1}{64}h^2 + \frac{1}{256}h^3 + \frac{25}{16384}h^4 + \frac{49}{65536}h^5 + \frac{441}{1048576}h^6 + \dots
$$

At $h = 1$ (where $b = 0$, a flat double line segment of total length $C(a, 0) = 4a$), the normalized circumference evaluates via Gauss's hypergeometric theorem to:

$$
{}_2F_1\!\left(-\frac{1}{2}, -\frac{1}{2}; 1; 1\right) = \frac{\Gamma(1)\,\Gamma(2)}{\Gamma(3/2)\,\Gamma(3/2)} = \frac{4a}{\pi a} = \frac{4}{\pi} \approx 1.2732395447351627
$$

---

## 2. Structural Proof of Ramanujan's Second Approximation

Historical reconstructions by Almkvist & Berndt (1988), Lorentzen (Jacobsen) & Waadeland (1992), and Villarino (2005/2008) noted that Ramanujan's perimeter approximations are deeply connected to continued fraction expansions and quadratic Padé approximants of Ivory's series. In `src/ramanujan/reconstruction.py`, we make the exact mechanism of the "4th-order miracle" explicit by converting $F(h) - 1 = {}_2F_1(-1/2, -1/2; 1; h) - 1$ into its canonical Stieltjes continued fraction (S-fraction):

$$
F(h) - 1 = \cfrac{a_1\,h}{1 - \cfrac{a_2\,h}{1 - \cfrac{a_3\,h}{1 - \cfrac{a_4\,h}{1 - \cfrac{a_5\,h}{1 - \ddots}}}}}
$$

Symbolic computation yields the exact rational sequence of partial coefficients:

$$
a_1 = \frac{1}{4}, \quad a_2 = \frac{1}{16}, \quad \mathbf{a_3 = \frac{3}{16}}, \quad \mathbf{a_4 = \frac{3}{16}}, \quad a_5 = \frac{11}{48}, \quad a_6 = \frac{115}{528}
$$

Because $a_3 = a_4 = 3/16$, replacing the continued fraction tail from depth 3 onwards by a constant 1-periodic continued fraction $T(h) = 1 - \frac{3h/16}{T(h)} = \frac{2 + \sqrt{4 - 3h}}{4}$ matches **both** depth 3 ($h^3$) and depth 4 ($h^4$) simultaneously and collapses into $1 + \frac{3h}{10 + \sqrt{4 - 3h}}$.

---

## 3. Inner-Radical Expansion and the Five 2D `HR` Approximations

Writing the exact normalized circumference in the generalized single-radical form

$$
F(h) = 1 + \frac{3h}{10 + \sqrt{4 - 3h - \delta(h)}} \iff \delta(h) = 4 - 3h - \left(\frac{3h}{F(h) - 1} - 10\right)^2
$$

and expanding $\delta(h)$ and its reciprocal $Q(h) = \frac{27\,h^4}{\delta(h)}$ symbolically around $h = 0$ reveals:

$$
Q(h) = \frac{27\,h^4}{\delta(h)} = 6144 - 7040\,h + \frac{5144}{3}\,h^2 - \frac{413}{9}\,h^3 - \frac{1081}{108}\,h^4 + \mathcal{O}(h^5)
$$

For the formula to be **100% exact at the degenerate boundary $h = 1$** ($F(1) = 4/\pi$), the inner correction $\delta(1) = \delta_1$ must satisfy:

$$
\delta_1 = 1 - \left(\frac{3\pi}{4 - \pi} - 10\right)^2 \approx 0.04082049800044723, \qquad Q(1) = \frac{27}{\delta_1} \approx 661.4324009399441
$$

### Taxonomy of Discovered Hyper-Ramanujan (`HR`) Formulas

1. **`HR-7-Int` (Pure Integer Order-7 Rational-Radical):**

$$
\frac{C}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{81\,h^4}{18432 - 21120\,h + 5144\,h^2}}}
$$

> Pure integers; exact up to $\mathcal{O}(h^8)$ at $h \to 0$; maximum error **78.05 ppm** at $h=1$.

2. **`HR-6B` (Boundary-Exact Order-6 Radical):**

$$
\frac{C}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\,h^4}{6144 - 7040\,h + B_{\text{exact}}\,h^2}}}, \qquad B_{\text{exact}} = \frac{27}{\delta_1} + 896 \approx 1557.43240094
$$

3. **`HR-8B` (Boundary-Exact Order-8 Radical):**

$$
\frac{C}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\,h^4}{6144 - 7040\,h + \frac{5144}{3}\,h^2 - \frac{413}{9}\,h^3 + C_{\text{exact}}\,h^4}}}, \qquad C_{\text{exact}} = \frac{27}{\delta_1} - \frac{6955}{9} \approx -111.34537684
$$

4. **`HR-8P` ★ (Singular-Pole Padé $-13/14$, Order-8 + Exact Boundary — Pareto Champion):**

$$
\frac{C}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27\left(1 - \frac{13}{14}h\right)h^4}{6144 - \frac{89216}{7}\,h + \frac{173288}{21}\,h^2 - \frac{103199}{63}\,h^3 + p_4\,h^4}}}, \qquad p_4 = \frac{27}{14\,\delta_1} - \frac{793}{63} \approx 34.65786991
$$

> Exact up to $\mathcal{O}(h^9)$ at $h \to 0$, exact at $h=1$, full 16-digit `float64` precision for $h \le 0.1$, and global maximum error **0.63 ppm** across $h \in [0, 1]$.

5. **`HR-6M` (Chebyshev-Minimax Singular-Padé $[3/1]$ + Exact Boundary):**

$$
\frac{C}{\pi(a+b)} \approx 1 + \frac{3h}{10 + \sqrt{4 - 3h - \dfrac{27(1 + q_1 h)h^4}{6144 + p_1 h + p_2 h^2 + p_3 h^3}}}
$$

> With $q_1 = -0.94542891$, $p_1 = -12848.71522304$, $p_2 = 8316.23093040$, and $p_3 = -1575.42062028$. Exact up to $\mathcal{O}(h^7)$ at $h \to 0$, exact at $h=1$, and global maximum error of **0.30 ppm** across $h \in [0, 1]$.

---

## 4. Pareto-Efficiency & Accuracy Benchmark

| Approximation Formula | Taylor Order | `sqrt` / `pow` | Op Cost | Digits (`h=0.1`) | Min Digits (`h <= 0.5`) | Min Digits (`h in [0,1]`) | Max Global Error (`ppm`) | **Pareto Efficiency Score** |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| **Kepler (1609)** | `h⁰` | `1` / `0` | `13.0` | 1.30 | 0.61 | 0.22 | `607300.92` | `0.214` |
| **Peano/Sipos (1792)** | `h¹` | `1` / `0` | `13.0` | 3.31 | 1.82 | 0.75 | `178097.25` | `0.582` |
| **Ramanujan I (1914)** | `h²` | `1` / `0` | `10.0` | 5.69 | 3.51 | 2.38 | `4155.03` | `1.279` |
| **Hudson [2/1] (1917)** | `h³` | `0` / `0` | `9.0` | 7.23 | 4.28 | 2.72 | `1889.83` | `1.665` |
| **Jacobsen-Waadeland [2/2] (1985)** | `O(h⁴)` | `0` / `0` | `13.0` | 8.85 | 5.12 | 3.05 | `888.05` | `1.662` |
| **Cantrell (2004, `+c·h¹²`)** | `O(h⁴)` | `1` / `0` | `23.5` | 9.58 | 5.79 | 4.86 | `13.74` | `1.449` |
| **Ramanujan II (1914)** | `O(h⁴)` | `1` / `0` | `15.0` | 9.58 | 5.76 | 3.40 | `402.34` | `1.698` |
| **Oertel (2026) [`HR-6B`]** | `O(h⁶)` | `1` / `0` | `28.5` | 13.13 | 7.73 | 4.93 | `11.82` | `1.695` |
| **Oertel (2026) [`HR-6M`]** | `O(h⁶)` | `1` / `0` | `37.5` | 13.61 | 8.34 | **6.52** | **`0.30`** | `1.618` |
| **Oertel (2026) [`HR-8B`]** | `O(h⁸)` | `1` / `0` | `32.5` | 15.31 | 8.48 | 5.04 | `9.11` | `1.786` |
| **Oertel (2026) [`HR-7-Int`]** | `O(h⁷)` | `1` / `0` | `26.0` | 14.64 | 8.40 | 4.11 | `78.05` | `1.885` |
| **Oertel (2026) [`HR-8P` ★]** | **`O(h⁸)`** | **`1` / `0`** | `34.5` | **16.00** | **9.81** | **6.20** | **`0.63`** | **`1.907` (Rank 1)** |

---

## 5. Generalization to 3D Triaxial Ellipsoid Surface Area (`HR-3D-10P`)

For a general triaxial ellipsoid with semi-axes $a \ge b \ge c \ge 0$, Carlson's theorem gives $S(a,b,c) = 4\pi\,R_G(X^2, Y^2, Z^2)$ with $X = ab$, $Y = ac$, $Z = bc$, $M = (X+Y+Z)/3$. In trace-free invariants $H = (x^2+y^2+z^2)/6 \in [0, 1]$ and $J = xyz \in [-1/4, 2]$, the exact 10th-order spherical series is:

$$
T_{10}(H, J) = 1 + \frac{3}{5}H - \frac{3}{35}J + \frac{3}{35}H^2 - \frac{6}{385}HJ + \frac{171}{5005}H^3 - \frac{3}{715}J^2 - \frac{3}{1001}H^2 J - \frac{57}{17017}H J^2 + \frac{27}{1547}H^4 - \frac{15}{46189}J^3 + \frac{36}{46189}H^3 J - \frac{6}{2717}H^2 J^2 + \frac{27}{2717}H^5
$$

On the entire degenerate 3D boundary $Z = bc \to 0$ (where $1 - 3H + J = 0$), $S(a,b,c)$ reduces identically to $\frac{\pi}{2}\,C_{\text{2D}}(X, Y)$. Using Spherical Ray Projection $s = 1 - Z/M \in [0, 1]$, $(H_b, J_b) = (H/s^2, J/s^3)$, and $h_{\text{ray}} = (4H_b - 1)/3 = ((X-Y)/(X+Y-2Z))^2 \in [0, 1]$ gives **`HR-3D-10P`**:

$$
S(a,b,c) \approx 4\pi\,M \left( T_{10}(H, J) + \frac{s^{11}}{9 - 8s}\left[\frac{3\pi}{8}\,F_{\text{2D, HR-8P}}(h_{\text{ray}}) - T_{10}(H_b, J_b)\right] \right)
$$
