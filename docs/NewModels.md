# New Models for VLEO Uncertainty Sources

This note summarizes (1) **existing** mathematical modeling approaches commonly used for each major Very Low Earth Orbit (VLEO) uncertainty source and (2) **new / stronger** modeling ideas intended to better capture **parametric**, **stochastic**, and **model-form** uncertainty in a way that is directly usable in **POD**, **orbit control**, **mission operations**, and **attitude control**.

---

## 0) Common backbone (applies to every source)

State-space form:
\[
\dot{\mathbf x}= \mathbf f(\mathbf x,\mathbf u,t;\boldsymbol\theta) + \mathbf w(t),
\qquad
\mathbf y_k=\mathbf h(\mathbf x(t_k),t_k;\boldsymbol\theta)+\mathbf v_k
\]
- \(\boldsymbol\theta\): uncertain parameters (static or slowly varying).
- \(\mathbf w(t)\): process noise (unmodeled accel/torque).
- \(\mathbf v_k\): measurement noise.

**Physics + discrepancy split (recommended for weakly-modeled sources):**
\[
\mathbf a(t)=\mathbf a_{\text{phys}}(t;\boldsymbol\theta)+\Delta\mathbf a(\mathbf s(t))
\]
where \(\mathbf s(t)\) are regime variables (e.g., altitude, latitude, local solar time, \(F_{10.7}\), \(A_p\), attitude angles, wall temperature).

**Two reusable modeling patterns:**

1) **Multiplicative log-error** for positive quantities
\[
q(t)=q_0(t)\exp(\delta_q(t)),\quad \delta_q \ \text{Gauss–Markov or GP}
\]

2) **Bounded parameters** via logistic transform
\[
p(t)=\frac{1}{1+e^{-z(t)}},\quad z(t)\ \text{Gauss–Markov}
\]

---

## 1) Thermospheric density uncertainty

### Existing modeling (typical)
- Use a semi-empirical model (NRLMSISE-00, JB2008, DTM family) and model residual as a **multiplicative bias/scale**:
\[
\rho(t)=\rho_{\text{model}}(t)\,\exp(\delta_\rho(t))
\]
- \(\delta_\rho(t)\) modeled as:
  - piecewise constant / piecewise linear per arc, or
  - 1st-order Gauss–Markov:
\[
d\delta_\rho = -\kappa_\rho \delta_\rho\,dt + \sigma_\rho\,dW_t
\]
- Assimilative correction fields (HASDM-like): global density corrections updated frequently from calibration satellites.

### New / stronger modeling ideas
1) **Spatio-temporal low-rank GP / basis for log-density error**
\[
\delta_\rho(\mathbf r,t)=\sum_{i=1}^p \beta_i(t)\,\phi_i(\mathbf r),
\qquad d\boldsymbol\beta = \mathbf A\boldsymbol\beta\,dt+\mathbf B\,d\mathbf W
\]
2) **Driver-aware stochastic density error**
\[
d\delta_\rho = -\kappa(\mathbf e)\delta_\rho\,dt + \sigma(\mathbf e)\,dW,\quad
\mathbf e=[F_{10.7},A_p,\text{LT},\text{lat},h]
\]
3) **Bayesian model averaging over multiple density models**
\[
\rho = \sum_m \pi_m \rho_m,\quad \pi_m\ge 0,\ \sum_m \pi_m=1
\]
with \(\pi_m\) adapted from residual likelihood.

---

## 2) Thermospheric winds uncertainty

### Existing modeling
- Use empirical wind climatology (e.g., HWM14):
\[
\mathbf V_{\text{rel}}=\mathbf V_{\text{inertial}}-\mathbf V_{\text{atm}}(\mathbf r,t)
\]
- OD practice often absorbs wind errors using empirical cross-track accelerations or inflated process noise.

### New / stronger modeling ideas
1) **Estimate wind as a stochastic low-rank field**
\[
\mathbf V_{\text{atm}}=\mathbf V_{\text{HWM}}+\sum_{i=1}^p \gamma_i(t)\,\boldsymbol\psi_i(\mathbf r)
\]
with \(\gamma_i(t)\) Gauss–Markov states.

2) **Activity-dependent switching (quiet vs disturbed)**
- Storm mode increases wind-error variance and shortens correlation time.

---

## 3) Composition uncertainty (O, N\(_2\), He, ...)

### Existing modeling
- Composition from thermosphere models influences mean molecular mass, partial densities, and energy exchange.
- Often **not estimated directly**; composition uncertainty is lumped into \(\delta_\rho\) and/or \(C_D\).

### New / stronger modeling ideas
1) **Simplex-respecting stochastic model for composition fractions**
Let \(\mathbf x_c \in \Delta^{K-1}\) (fractions sum to 1). Use logistic-normal:
\[
\mathbf z\in\mathbb R^{K-1},\quad \mathbf x_c=\text{softmax}([\mathbf z,0])
\]
\[
d\mathbf z = -\mathbf K(\mathbf z-\boldsymbol\mu)\,dt+\mathbf\Sigma^{1/2}\,d\mathbf W
\]
2) **Joint density–composition estimation**
Estimate \(\delta_\rho\) and a reduced composition mode (e.g., \(O/N_2\)).

---

## 4) Gas–surface interaction (GSI) / drag coefficient uncertainty

### Existing modeling
**(a) Reduced-coefficient approach**
\[
\mathbf a_D = -\frac{1}{2}\frac{\rho}{m} C_D A\,\|\mathbf V_{\text{rel}}\|\,\mathbf V_{\text{rel}}
\]
with \(C_D\) or ballistic coefficient \(B=C_DA/m\) estimated per arc (piecewise constant) or random walk.

**(b) Physical / semi-physical GSI models**
Maxwell/Sentman/DRIA/CLL-like kernels; DSMC/TPMC for complex geometry; surrogate models for fast evaluation.

### New / stronger modeling ideas
1) **Kernel-parameter state model (bounded Gauss–Markov)**
For CLL-like \(\alpha_n,\alpha_t\in(0,1)\):
\[
\alpha=\text{logistic}(z),\quad dz = -\kappa(z-\mu)\,dt+\sigma\,dW
\]
Panel-wise: \(z_j=z_{\text{common}}+\tilde z_j\) to capture global aging + local variation.

2) **Model-choice uncertainty as a mixture of kernels**
\[
K = \lambda K_{\text{CLL}}+(1-\lambda)K_{\text{Maxwell}},\quad \lambda\in[0,1]
\]
Estimate \(\lambda\) as a slow random-walk state.

3) **AO-fluence-driven drift**
\[
\Phi(t)=\int n_{AO}(t)\|\mathbf V_{\text{rel}}\|\,dt,\quad
\alpha(t)=g(\Phi(t))+\epsilon(t)
\]
Links GSI uncertainty directly to atomic oxygen exposure.

4) **Explicit force/torque discrepancy**
\[
\mathbf a_{\text{aero}} = \mathbf a_{\text{phys}}(\boldsymbol\phi) + \Delta\mathbf a(\mathbf s)
\]
with \(\Delta\mathbf a\) as low-rank basis or GP:
\[
\Delta\mathbf a(\mathbf s)=\mathbf B(\mathbf s)\boldsymbol\beta,\quad \boldsymbol\beta\sim\mathcal N(0,\mathbf\Sigma_\beta)
\]

---

## 5) Rarefied-regime / Knudsen effects

### Existing modeling
- Assume free molecular flow above some altitude; use analytic FMF theory for simple shapes; DSMC/TPMC otherwise.
- Transitional regimes sometimes handled with bridging relations.

### New / stronger modeling ideas
1) **Latent regime weight with uncertain transition altitude**
\[
C = w(Kn;\,Kn_0)\,C_{\text{FMF}} + (1-w)\,C_{\text{cont}}
\]
where \(Kn_0\) is uncertain/estimated; \(w(\cdot)\) can be logistic.

2) **Multi-fidelity discrepancy learning**
Learn \(\Delta C(Kn,\alpha,\beta,T_w,\dots)\) from sparse DSMC truth using GP / polynomial chaos.

---

## 6) Space-weather drivers (F10.7, Kp/Ap) & storm-time disturbances

### Existing modeling
- Indices are treated as deterministic inputs to density/wind models.
- Assimilation systems update corrections frequently and use short-term forecasting.

### New / stronger modeling ideas
1) **Treat drivers as uncertain exogenous states**
\[
d\mathbf e = \mathbf a(\mathbf e)\,dt + \mathbf B\,d\mathbf W
\]
Propagate driver uncertainty into density/wind predictions.

2) **Change-point / switching storm model**
Hidden mode \(m(t)\in\{\text{quiet},\text{storm}\}\) switches \(\sigma_\rho,\sigma_V\) and bias dynamics.

---

## 7) SRP and Earth albedo/IR uncertainties

### Existing modeling
- Cannonball or box-wing SRP; estimate optical coefficients.
- Hybrid physical + empirical periodic terms (e.g., ECOM-style).

### New / stronger modeling ideas
1) **Optical coefficient aging as bounded stochastic states**
\[
\rho_{\text{opt}}=\text{logistic}(z),\quad dz=-\kappa(z-\mu)\,dt+\sigma\,dW
\]
2) **Sun-angle / eclipse structured discrepancy**
\[
\Delta \mathbf a_{\text{SRP}}=\mathbf B(\epsilon,\text{eclipse})\boldsymbol\beta
\]
3) **Joint SRP + thermal recoil modeling** (if thermal telemetry exists)

---

## 8) Gravity field, tides/loading, and third-body ephemerides

### Existing modeling
- Gravity: spherical-harmonic geopotential models with published coefficient covariances / degree variances.
- Tides/loading: IERS Conventions (solid Earth tides, pole tides, ocean loading, etc.).
- Third-body: JPL DE ephemerides.

### New / stronger modeling ideas
1) **Compress gravity-field uncertainty into colored acceleration noise**
Match low-order colored process statistics to degree-variance-induced accel error.

2) **Time-variable gravity residual as seasonal + colored noise**
Useful for long arcs or when modeling time-variable gravity imperfectly.

3) **Selective ephemeris uncertainty propagation**
Include only for high-precision covariance analysis; otherwise negligible in many VLEO contexts.

---

## 9) Magnetic field and residual dipole torques

### Existing modeling
Attitude torque:
\[
\boldsymbol\tau_m = \mathbf m_{\text{res}} \times \mathbf B(\mathbf r,t)
\]
- \(\mathbf B\) from IGRF-type models.
- \(\mathbf m_{\text{res}}\) treated as constant unknown or slowly varying (sometimes calibrated).

### New / stronger modeling ideas
1) **Augment \(\mathbf m_{\text{res}}\) as stochastic state**
\[
d\mathbf m = -\kappa(\mathbf m-\bar{\mathbf m})\,dt+\Sigma^{1/2}d\mathbf W
\]
(optionally temperature-dependent \(\bar{\mathbf m}(T)\)).

2) **Explicit field discrepancy**
\[
\mathbf B=\mathbf B_{\text{IGRF}}+\delta \mathbf B
\]
Use Gaussian or heavy-tail models for \(\delta\mathbf B\) if residuals are non-Gaussian.

---

## 10) Atomic oxygen (AO) erosion impacts

### Existing modeling
- Often handled qualitatively or via material degradation/erosion models.
- In aero: AO changes roughness and accommodation, thus modifying \(C_D\)/GSI parameters.

### New / stronger modeling ideas
1) **AO fluence as mission state driving parameter drift**
\[
\Phi(t)=\int \rho_{AO}\|\mathbf V_{\text{rel}}\|\,dt
\]
Let \(\alpha(\Phi)\), reflectivity(\(\Phi\)), roughness(\(\Phi\)) evolve with uncertainty bands.

2) **Switched “surface condition” model**
\[
s(t)\in\{\text{fresh},\text{aged},\text{contaminated}\}
\]
Mode affects GSI and optical properties; transitions depend on events (maneuvers/outgassing) and environment.

---

## 11) Propulsion/thrust and mass properties

### Existing modeling
Scale, bias, misalignment:
\[
\mathbf F_T = \eta\,\mathbf F_{\text{cmd}} + \mathbf b + \mathbf n(t),\qquad
\mathbf T_T = \mathbf r_T\times \mathbf F_T
\]
- \(\eta\): constant or slow drift.
- \(\mathbf b\): bias (often random walk).
- Misalignment modeled as small-angle uncertainty.
- Mass \(m(t)\) propagated from prop usage with uncertainty.

### New / stronger modeling ideas
1) **Telemetry-informed stochastic actuator model**
Make \(\eta(T,p)\) depend on temperature/pressure states and estimate online.

2) **Impulse-bit heavy-tail model**
Model delivered impulse with mixture/heavy-tail to capture anomalies and valve stiction.

---

## 12) Measurement/sensor errors (GNSS, SLR, attitude sensors/IMU, clocks)

### Existing modeling
- GNSS: white noise + elevation/CN0-dependent variance; biases (multipath, phase ambiguity).
- SLR: single-shot mm–cm; normal points often mm-level depending on station/product.
- IMU: ARW, bias instability, RRW modeled as random walk / Gauss–Markov; identified via Allan variance.
- Clock: random walk / GM in many OD/PPP formulations.

### New / stronger modeling ideas
1) **Robust measurement noise (mixtures/heavy tails)**
\[
\mathbf v_k \sim (1-\epsilon)\mathcal N(0,R)+\epsilon\,t_\nu(0,R')
\]
2) **Time-correlated measurement noise**
ARMA / multi-GM processes for clocks and IMU biases.

3) **Adaptive covariance from diagnostics**
Update \(R(t)\) using CN0, eclipse flags, thermal states, tracking geometry.

---

## 13) Operational dropouts/latency and software/human factors

### Existing modeling
- Often handled procedurally (gaps) rather than probabilistically.
- OD: missing observations; ops: contingency margins.

### New / stronger modeling ideas
1) **Dropout model (Bernoulli / Markov)**
\[
\Pr(\text{obs available at }k)=p_k,\quad p_k=\sigma(\mathbf q^\top \mathbf z_k)
\]
where \(\mathbf z_k\) contains link margin, antenna geometry, eclipse, ops mode.

2) **Random delay model for latency**
\[
\tau_k \sim \text{LogNormal}(\mu_\tau,\sigma_\tau)
\]
Handle via out-of-sequence measurement filtering.

3) **Hybrid-system mode state for ops anomalies**
\[
m(t)\in\{\text{nominal},\text{safe},\text{degraded sensors},\text{ops anomaly}\}
\]
Mode switches process/measurement models (effective for mission control realism).

---

## Practical workflow to build & test “new models”

1) Pick one source (high payoff: **GSI**, **storm-time density/winds**, **aero torque / CP drift**).
2) Choose representation:
   - augmented parameter state (GM / RW),
   - structured discrepancy \(\Delta\mathbf a(\mathbf s)\),
   - hybrid switching model.
3) Write the discrete-time version for your estimator/controller:
   - state augmentation,
   - transition Jacobian \(F\),
   - process covariance \(Q\),
   - measurement sensitivities.
4) Define identifiability hooks:
   - which measurements constrain which uncertainty (GNSS/SLR residual signatures, accelerometer, attitude rates, etc.).

---

## Reference pointers (starting points)

- Density rescaling / calibration in OD; thermosphere model performance benchmarking:
  - https://www.swsc-journal.org/articles/swsc/full_html/2012/01/swsc120006/swsc120006.html
- HASDM / assimilation concept (density corrections):
  - https://www.sciencedirect.com/science/article/pii/S0273117705002048
- HWM14 wind model description:
  - https://ccmc.gsfc.nasa.gov/models/HWM14~2014/
- Drag/GSI sensitivity in free molecular flow contexts:
  - https://www.sciencedirect.com/science/article/pii/S027311771400413X
- Comparing physical drag coefficients; DRIA/CLL-style ideas:
  - https://laro.lanl.gov/esploro/outputs/journalArticle/Comparing-Physical-Drag-Coefficients-Computed-Using/9916369938803761
- SRP “adjustable box-wing” approach:
  - https://portal.fis.tum.de/en/publications/adjustable-box-wing-model-for-solar-radiation-pressure-impacting-/
- IERS Conventions (tides/loading standards):
  - https://www.iers.org/IERS/EN/Publications/TechnicalNotes/tn36
- JPL DE ephemerides:
  - https://ssd.jpl.nasa.gov/doc/de430_de431.html
- IGRF overview:
  - https://geomag.bgs.ac.uk/research/modelling/IGRF
- GNSS error characterization example:
  - https://www.mdpi.com/1424-8220/20/14/4046
- IMU error modeling and Allan-variance context:
  - https://jeas.springeropen.com/articles/10.1186/s44147-024-00520-9