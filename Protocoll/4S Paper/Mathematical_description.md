# Mathematical and Quantitative Descriptions of Uncertainty Sources in VLEO

## Executive summary

Very Low Earth Orbit (VLEO, here taken as 150–450 km) is the regime where **the dominant uncertainty in orbit dynamics is usually aerodynamic**, because neutral density, composition, and winds vary strongly with solar/geomagnetic forcing and with local time, and the mapping from environment to force/torque depends on uncertain rarefied gas–surface interaction (GSI) physics. citeturn11view0turn0search0turn8search35

Across the literature, uncertainty for POD/estimation and propagation is typically represented using a small set of **stochastic primitives**:

- **Additive white noise** (often Gaussian) for measurement noise in GNSS/SLR/attitude sensors. citeturn4search24turn3search16turn3search29  
- **Additive constant bias + random walk** for slowly drifting sensor/clock biases and some model parameters (e.g., accelerometer bias, gyro bias, receiver clock). citeturn6search24turn3search11turn3search17  
- **First‑order Gauss–Markov (exponentially correlated) processes** for unmodelled accelerations (dynamic model compensation, DMC) and for time-correlated drag/radiation-pressure residuals used to obtain covariance realism in OD. citeturn6search13turn7search29turn8search20  
- **Multiplicative scale factors** (often modelled in log space) for density, drag coefficient, SRP coefficient, thrust magnitude, and radiation-pressure scaling. citeturn12search18turn7search0turn7search25  
- **Non‑Gaussian / heavy‑tailed mixtures** for outliers (tracking multipath spikes, cycle slips, SLR detector artefacts), operational anomalies, and human/software-induced errors, often handled with robust estimation or mixture noise models rather than pure Gaussian assumptions. citeturn3search20turn4search16turn7search18

Representative quantitative ranges that are widely used for VLEO modelling include: neutral mass density spanning roughly **10⁻¹⁰–10⁻¹³ kg m⁻³** across 200–450 km depending on solar activity and diurnal state; “quiet” thermosphere density modelling uncertainty on the order of **~10–15%**; and storm-time deviations that can be much larger. citeturn11view0turn12search18turn6search12

## Notation and stochastic-process building blocks

### Dynamics and observation equations with explicit uncertainty terms

Orbit/attitude estimation commonly uses a stochastic state-space model of the form

\[
\dot{\mathbf{x}}(t)=\mathbf{f}\!\left(\mathbf{x},t;\boldsymbol{\theta}\right)+\mathbf{G}(t)\,\mathbf{w}(t),
\qquad
\mathbf{y}_k=\mathbf{h}\!\left(\mathbf{x}(t_k),t_k;\boldsymbol{\eta}\right)+\mathbf{b}(t_k)+\boldsymbol{\nu}_k,
\]

where \(\boldsymbol{\theta}\) are uncertain force/torque parameters (density scale factors, \(C_D\), \(C_R\), thrust scale, etc.), \(\mathbf{w}(t)\) captures unmodelled dynamics (often coloured), and \(\mathbf{b},\boldsymbol{\nu}_k\) are measurement biases and noise. DMC/SNC formulations explicitly motivate modelling \(\mathbf{w}(t)\) with Gauss–Markov structure to obtain realistic covariances. citeturn6search13turn7search29turn6search16

A standard decomposition for non-gravitational accelerations is

\[
\ddot{\mathbf{r}}= \mathbf{a}_{grav}(\mathbf{r},t)
+\mathbf{a}_{drag}(\rho, \mathbf{u},C_D, \mathbf{q})
+\mathbf{a}_{RP}(C_R, \text{optics}, \mathbf{q})
+\mathbf{a}_{th}(T,\mathbf{u}_T)
+\mathbf{a}_{emp}(t),
\]

with the key VLEO sensitivity being \(\mathbf{a}_{drag}\) and its parameter/environment uncertainty. citeturn11view0turn8search35turn7search25

### Standard stochastic models encountered in OD/AOCS literature

**White Gaussian noise (WGN):**
\[
\nu_k \sim \mathcal{N}(0, R_k)
\]
Used for many measurement models, sometimes inflated or elevation-dependent for GNSS. citeturn4search24turn4search36turn3search11

**Random walk (RW):**
\[
b_{k+1}=b_k+\eta_k,\quad \eta_k\sim \mathcal{N}(0, q\,\Delta t)
\]
Used for slowly drifting biases and certain coefficients treated as “nearly constant but drifting”. citeturn6search24turn7search14turn3search17

**First-order Gauss–Markov / Ornstein–Uhlenbeck (GM(1)):**
\[
\dot{z}(t)= -\frac{1}{\tau}z(t)+\sigma\sqrt{\frac{2}{\tau}}\,w(t)
\]
Stationary variance \(\mathrm{Var}[z]=\sigma^2\), correlation \(C(\Delta t)=\sigma^2 e^{-|\Delta t|/\tau}\). This is the canonical model for exponentially correlated unmodelled accelerations in dynamic model compensation. citeturn6search13turn7search29turn6search16

**Discrete AR(1) form of GM(1):**
\[
z_{k+1}=\phi z_k+\eta_k,\quad \phi=e^{-\Delta t/\tau},\quad
\eta_k\sim\mathcal{N}\!\left(0,\sigma^2(1-\phi^2)\right)
\]
Used directly in discrete Kalman filtering. citeturn6search16turn6search24turn8search20

**Multiplicative/lognormal scale factor models** (common for density and coefficients):
\[
\rho(t)=\rho_{model}(t)\,\exp(\delta_\rho(t)), \quad \delta_\rho \sim \text{GM(1) or RW}
\]
This avoids negative densities and naturally encodes “percentage” errors. citeturn12search18turn8search1turn12search19

### Mermaid diagram for common temporal correlation structures

```mermaid
flowchart TD
  W[white noise w(t)] --> RW[random walk: b_{k+1}=b_k+η]
  W --> GM[first-order Gauss–Markov: ż=-(1/τ)z+σ*sqrt(2/τ)w]
  GM --> AR1[discrete AR(1): z_{k+1}=φ z_k+η]
  RW --> IRW[integrated random walk (clock bias)]
  GM --> GM2[2nd-order Gauss–Markov (oscillatory / coloured)]
```

Time-correlated error modelling is explicitly emphasised in OD/process-noise estimation work as central to covariance realism in LEO operations. citeturn6search13turn7search29turn12search32

## Mathematical and quantitative models for each VLEO uncertainty source

The entries below focus strictly on **uncertainty representation** (stochastic/statistical model, parameter ranges, correlations, and how the term enters the equations). Where the literature reports ranges outside 150–450 km, values are interpreted for the VLEO band when feasible.

**Neutral thermospheric density \(\rho\)**  
*Uncertainty model:* a common OD representation is a **multiplicative density scale factor** relative to an empirical model (e.g., MSIS/DTM), modelled as a constant over an arc, as a RW, or as a GM(1) process:
\[
\rho(t)=\rho_{mod}(t)\,s_\rho(t),\quad s_\rho(t)=\exp(\delta_\rho(t)),\quad \delta_\rho \sim \text{GM(1) or RW}.
\]
Density corrections (“scale factors”) are reported and compared across techniques (e.g., SLR vs accelerometer), often with multi-hour cadence (e.g., 12 h) in published products. citeturn8search1turn12search18turn12search17  
*Representative ranges (150–450 km):* a classic NASA environment guideline figure shows order-of-magnitude densities spanning, for example, \(\sim10^{-9}\) kg m⁻³ near 200 km and \(\sim10^{-12}\)–\(10^{-13}\) kg m⁻³ near 400–450 km depending on solar activity and diurnal state (plot axis in g cm⁻³ converted to SI). citeturn11view0 Density-model uncertainty under **quiet conditions** is often stated around **10–15%** (order-of-magnitude) for operational thermosphere modelling, while much larger deviations occur during storms. citeturn12search18turn6search12  
*Correlation scales:* OD literature increasingly stresses that density error is **time-correlated** (not white); operationally relevant correlation “half-life” / \(\tau\) is treated as a tuning/estimated parameter in stochastic drag models. citeturn12search26turn12search23turn12search16 Spatial correlation is anisotropic; practical along-track aggregation scales on the order of **~10,000 km** appear in uncertainty analyses of thermospheric mass density observations derived from tracking/accelerometers (reflecting both physics and retrieval noise structure). citeturn12search3turn8search20  
*Entry into POD/propagation:* via drag,
\[
\mathbf{a}_D=-\tfrac12 \rho C_D \frac{A_\perp}{m}\,|\mathbf{v}_{rel}|\mathbf{v}_{rel}.
\]
An uncertain \(\rho\) is thus a multiplicative uncertainty on \(\mathbf{a}_D\). citeturn12search19turn8search35

**Thermospheric winds \(\mathbf{u}\)**  
*Uncertainty model:* winds typically enter as uncertainty in relative velocity \(\mathbf{v}_{rel}=\mathbf{v}-\mathbf{u}\). OD/FD often represents wind error as (i) additive WGN on wind components, (ii) a correlated GM(1) on wind components, or (iii) as part of the effective drag-acceleration residual absorbed by empirical accelerations. citeturn0search1turn6search13turn0search32  
*Representative ranges:* empirical wind models (e.g., HWM14) provide climatological winds; discrepancies between HWM14 and measurements can be substantial seasonally/regionally, consistent with wind variability at thermospheric heights. citeturn0search1turn0search24 NASA’s orbital environment guideline notes thermospheric circulation with “maximum velocities … typically hundreds of metres per second.” citeturn9view0  
*Correlation scales:* winds have strong diurnal/tidal signatures in empirical models (migrating tides explicitly noted in HWM14 formulation); correlation times of hours are commonly implied for “quiet-time” structures, while shorter scales (minutes) arise from waves and disturbances. citeturn0search1turn9view0  
*Entry into POD/propagation:* through \(|\mathbf{v}-\mathbf{u}|\) and direction in \(\mathbf{a}_D\), and through aerodynamic lift/torque models that depend on flow incidence. citeturn0search32turn8search35

**Composition variability (e.g., \(\Sigma O/N_2\), mean molecular mass \(\bar{m}\))**  
*Uncertainty model:* composition enters density models and GSI physics; in OD it is often represented indirectly as part of a density scale factor, or explicitly by augmenting the state with composition-related parameters and treating them as GM(1)/RW processes. Composition ratios can be modelled as lognormal random fields or as deviations from a climatology with correlated noise. citeturn0search0turn5search4turn12search19  
*Representative ranges:* column \(\Sigma O/N_2\) ratios reported from FUV-based products commonly vary on regional and storm-time scales; published studies report typical mid-latitude values on the order of ~0.4–0.8 depending on conditions and interpretation. citeturn5search28turn5search4  
*Correlation scales:* composition responds to geomagnetic forcing and transport; storm-time composition depletion patterns evolve over hours to days and can extend from high to mid/low latitudes. citeturn5search0turn5search35  
*Entry into POD/propagation:* via \(\rho(\bar{m},T)\) in thermosphere models and via GSI coefficients \(C_D(\text{species},T_w)\); in joint estimation, composition uncertainty is strongly correlated with \(C_D\) and \(\rho\). citeturn8search35turn12search19

**Gas–surface interaction parameters (drag coefficient \(C_D\), accommodation, scattering kernels)**  
*Uncertainty model:* literature distinguishes fixed/fitted/physical approaches; uncertainty is represented as (i) a stochastic parameter \(C_D(t)\) estimated with RW/GM(1), (ii) uncertain GSI parameters (e.g., accommodation coefficients \(\alpha\)) propagated via a physical \(C_D(\alpha,\dots)\), or (iii) model-form uncertainty across scattering kernels (Maxwell, CLL, etc.) treated as discrete alternatives or mixtures. citeturn6search19turn8search35turn8search19  
A common fitted-parameter stochastic model is:
\[
C_{D,k+1}=C_{D,k}+\eta_k,\quad \eta_k\sim\mathcal{N}(0,q_{C_D}\Delta t)
\quad \text{or}\quad
\dot{C}_D=-(1/\tau) (C_D-\bar{C}_D)+\sigma w.
\]
Random-walk drag/SRP coefficients are explicitly documented in real-time POD configurations for LEO. citeturn7search14turn7search5  
*Representative ranges:* published operational tuning for a LEO mission reported exploring \(C_D\) roughly 1.6–2.5. citeturn7search5 VLEO-focused rarefied modelling studies report converged \(C_D\) values around ~1.2–1.6 for CubeSat-like orientations at higher altitudes, with sensitivity to specularity/accommodation possibly pushing \(C_D\) above ~2 in some cases. citeturn8search2turn8search35  
*Correlation scales:* \(C_D\) temporal correlation is driven by (a) attitude/orientation, (b) thermospheric composition, and (c) surface state evolution; stochastic “half-life” assumptions (hours–days) are used in modern uncertainty realism work. citeturn12search26turn6search19  
*Entry into POD/propagation:* multiplicatively in drag and also in torque models; for a macro‑model, force/torque coefficients depend on scattering kernel parameters. citeturn8search26turn8search19

**Rarefied flow / Knudsen regime uncertainty (Kn = \(\lambda/L\))**  
*Uncertainty model:* treated primarily as **model-form uncertainty** because regime classification (continuum ↔ transitional ↔ free-molecular) determines which aerodynamic model is valid. A practical uncertainty representation is to treat \(Kn\) (through \(\lambda\), hence \(\rho\) and composition) as uncertain and propagate it into coefficient uncertainty: \(C_D = C_D(Kn,\alpha,\dots)\) with uncertain inputs. citeturn8search25turn8search33turn8search22  
*Representative ranges:* VLEO experimental/analysis literature uses Knudsen number explicitly to characterise rarity and connects it to mean free path and spacecraft scale. citeturn8search25turn8search33 A VLEO drag modelling study reports Knudsen numbers and corresponding \(C_D\) behaviour across 50–500 km, illustrating that “regime-aware” modelling is needed. citeturn8search2  
*Correlation scales:* regime variation is strongly altitude-dependent (scale height), with shorter-term modulation through \(\rho\) variability (minutes–hours). citeturn11view0turn12search18  
*Entry into POD/propagation:* via selecting/parameterising \(\mathbf{a}_{drag}\) and aerodynamic torque models; uncertainty can be expressed as switching uncertainty or as inflated coefficient covariance when \(Kn\) is near the transitional regime. citeturn8search19turn8search22

**Atomic oxygen (AO) flux and surface erosion**  
*Uncertainty model:* AO flux to a ram-facing surface is often modelled as
\[
\Phi_{AO}(t)=n_{AO}(t)\,|\mathbf{v}_{rel}(t)|,
\]
with \(n_{AO}\) derived from atmospheric models or measurements; uncertainty is thus often multiplicative (through \(n_{AO}\)) plus orientation uncertainty (through projected area and incidence). citeturn2search6turn2search3  
Erosion/oxidation is frequently modelled with yield \(Y\) and fluence \(F=\int \Phi\,dt\), giving thickness/mass loss \(\Delta h \propto YF\); uncertainty enters through \(Y\) (material-dependent) and \(F\). citeturn2search3turn2search6  
*Representative ranges:* NASA AO environment references give incident AO fluxes around **10¹⁴–10¹⁵ atoms cm⁻² s⁻¹** in typical LEO altitude ranges (hundreds of km). citeturn2search6turn2search3  
*Correlation scales:* AO density follows thermosphere variability (diurnal to storm-time). Mission-level AO fluence specification explicitly includes uncertainty components in space environment specification practice. citeturn2search25turn12search18  
*Entry into POD/propagation:* primarily indirect, via time evolution of surface roughness/chemistry affecting \(C_D\), accommodation, and optics (thus both drag and radiation pressure coefficients drift). citeturn8search35turn2search3

**Plasma interactions and charging**  
*Uncertainty model:* charging is often treated with environment-driven uncertain currents and potentials (spacecraft potential \(V_s\)) governed by current balance; uncertainty representations include bounded uncertain plasma parameters (density/temperature), stochastic “event” processes for auroral precipitation, and random-walk biases for surface potentials during eclipse transitions. citeturn9view2turn0search25  
*Representative ranges (VLEO):* in LEO, high plasma density tends to keep potentials relatively low compared to high-altitude regimes, but the handbook emphasises that energetic particle environments (notably auroral) can produce charging concerns above ~250 km and during disturbed conditions. citeturn9view0turn9view2  
*Correlation scales:* charging can change rapidly (seconds–minutes) with environment changes and sunlight/eclipse; disturbance-driven intervals are often modelled as piecewise regimes. citeturn9view2turn0search25  
*Entry into POD/propagation:* usually indirect (sensor outages/biases, attitude disturbances) rather than as a dominant translational acceleration term in standard OD; in specialised models, plasma drag/interaction may be incorporated as an additional uncertain acceleration. citeturn9view2turn0search25

**Aerodynamic forces and torques (stochastic components beyond deterministic drag)**  
*Uncertainty model:* OD frequently represents remaining mismodelled aerodynamic accelerations (and sometimes torque disturbances) using empirical accelerations with GM(1) dynamics:
\[
\dot{\mathbf{a}}_{emp}=-(1/\tau)\mathbf{a}_{emp}+\mathbf{w}(t),
\]
with \(\tau\) chosen/tuned (minutes to one orbit) and \(\mathbf{w}\) white. citeturn6search13turn7search29turn6search16  
Published uncertainty analyses for thermosphere-derived accelerations explicitly adopt correlation times that differ by direction, e.g., along‑track correlation times on the order of **one orbital period (~5400 s)** and shorter times (hundreds of seconds) in radial/cross-track for LEO aerodynamic acceleration modelling. citeturn8search20turn12search3  
*Representative magnitudes:* unmodelled residual accelerations in simulation-based OD studies are often assumed at \(\lesssim10^{-7}\,\mathrm{m/s^2}\) scale, but in VLEO the aerodynamic acceleration itself can be much larger; hence the “residual” magnitude depends strongly on modelling fidelity. citeturn7search2turn11view0  
*Correlation scales:* density/wind-driven “weather” yields coloured acceleration noise; modern OD realism work targets modelling density auto-correlation explicitly rather than “white process noise”. citeturn12search26turn12search23  
*Entry into POD/propagation:* additively as \(\mathbf{a}_{emp}\) in the equations of motion and as torque disturbances in attitude dynamics. citeturn6search13turn7search29

**Solar EUV/UV and space‑weather indices (F10.7, Kp, Ap) uncertainty**  
*Uncertainty model:* indices enter density/wind models as drivers; uncertainty is represented as (i) measurement/calibration uncertainty on observed indices and (ii) forecast uncertainty for predicted indices used in propagation. Forecast error is often treated as an additive error on the index, sometimes with time correlation (e.g., AR(1)). citeturn1search1turn1search9turn12search18  
*Representative ranges/units:* F10.7 is in solar flux units (sfu); the measurement and calibration process for the 10.7 cm flux is described in detail in the primary reference paper. citeturn1search1 Kp is a 3‑hourly planetary index derived from standardised observatory K indices; Ap is a daily linearised index derived from 3‑hourly ap values with units effectively tied to **2 nT** steps as documented in GFZ descriptions. citeturn1search18turn1search22  
*Correlation scales:* Kp/ap are inherently piecewise-constant over 3‑hour bins; F10.7 has strong day-to-day correlation and 27‑day solar rotation structure, with forecast verification products documenting lead-time dependence. citeturn1search9turn1search1  
*Entry into POD/propagation:* through density model inputs \(\rho_{mod}(F10.7,ap,\ldots)\), so index uncertainty becomes parametric uncertainty in \(\rho\) and thus in \(\mathbf{a}_D\). citeturn0search0turn12search18

**Geomagnetic storm effects (thermosphere disturbance response)**  
*Uncertainty model:* storms are often modelled as regime-switching increases in density and winds; a convenient stochastic representation is a **jump or piecewise-bias component** added to the log-density error:
\[
\delta_\rho(t)=\delta_{\rho,\text{quiet}}(t)+\sum_i J_i\,\mathbf{1}_{[t_i,t_i+\Delta_i]}(t),
\]
where \(J_i\) are storm-time log-density anomalies. citeturn6search12turn6search7  
*Representative magnitudes:* in an orbital environment guideline, simulation results for a substorm suggest a significant **global density response (~60%)** at 200–450 km under that scenario. citeturn9view0 Modern model assessment papers define storm conditions via ap thresholds (e.g., ap ≥ 80) and evaluate storm-time performance/uncertainty accordingly. citeturn6search12turn1search22  
*Correlation scales:* storm-time density disturbances persist for hours to days and are spatially global/hemispheric; hence both temporal and spatial correlation lengths increase relative to quiet-time fluctuations. citeturn6search12turn9view0  
*Entry into POD/propagation:* primarily through \(\rho(t)\) and \(\mathbf{u}(t)\) affecting \(\mathbf{a}_D\), plus indirect impacts through GNSS performance degradation under severe events (treated as increased measurement noise/outages). citeturn6search12turn12search36

**Gravity field (higher-order geopotential) coefficient uncertainty**  
*Uncertainty model:* geopotential is expanded in spherical harmonics with uncertain coefficients \(C_{\ell m},S_{\ell m}\) with published error statistics/covariance. A standard linearised uncertainty injection is:
\[
C_{\ell m}=C_{\ell m}^0+\delta C_{\ell m},\quad \mathbf{\delta a}_{grav}\approx
\sum_{\ell,m}\frac{\partial \mathbf{a}_{grav}}{\partial C_{\ell m}}\,\delta C_{\ell m}
+\frac{\partial \mathbf{a}_{grav}}{\partial S_{\ell m}}\,\delta S_{\ell m}.
\]
EGM2008 was developed with associated error covariance information and is widely used as a baseline. citeturn1search0turn1search12  
*Representative ranges:* coefficient standard errors depend strongly on \(\ell,m\) and geography; studies explicitly use EGM2008 coefficient standard errors to build global covariance functions and then calibrate/scale them. citeturn1search28turn1search0  
*Correlation scales:* static gravity coefficient errors are temporally constant (static model), but map into orbit errors with characteristic orbital-period signatures; spatially, they correspond to wavelength \(\sim 2\pi R_E/\ell\). citeturn1search0turn1search32  
*Entry into POD/propagation:* via \(\mathbf{a}_{grav}\) in the dynamics; gravity-field uncertainty is especially relevant for radial/in-track precision when drag is well estimated. citeturn1search0turn12search3

**Earth tides and loading (solid Earth tide, ocean tide, pole tide, atmospheric loading)**  
*Uncertainty model:* most OD uses deterministic IERS-convention models; uncertainty can be described as parametric error in Love/Shida numbers or tide model amplitudes/phases:
\[
\Delta \mathbf{r}_{tide}(t)=\Delta\mathbf{r}_{IERS}(t;\mathbf{p})+\mathbf{e}_{tide}(t),
\quad \mathbf{e}_{tide}\ \text{often treated as coloured residual}.
\]
The IERS Conventions provide the canonical modelling equations and parameters for consistent practice. citeturn0search10turn0search26turn0search14  
*Representative magnitudes:* tidal site displacements are at the cm to dm level (geodetic context) and must be corrected for mm–cm orbit/SLR consistency; remaining modelling errors contribute to centimetre-to-millimetre residual structures depending on technique and station modelling. citeturn3search1turn0search26  
*Correlation scales:* tidal errors are periodic at known constituent frequencies; loading effects have meteorological correlation times (hours–days). citeturn0search26turn3search1  
*Entry into POD/propagation:* via station coordinate corrections in observation models (GNSS/SLR) and via time-variable geopotential terms. citeturn0search14turn3search1

**Third‑body ephemeris uncertainty (Sun/Moon/planets positions)**  
*Uncertainty model:* third-body acceleration for body \(b\) is
\[
\mathbf{a}_b=\mu_b\left(\frac{\mathbf{r}_b-\mathbf{r}}{|\mathbf{r}_b-\mathbf{r}|^3}-\frac{\mathbf{r}_b}{|\mathbf{r}_b|^3}\right).
\]
Uncertainty enters through ephemeris position errors \(\delta \mathbf{r}_b\) and \(\delta(\mu_b)\); ephemeris products provide accuracy context and sometimes covariance/uncertainty discussion. citeturn5search1turn5search21  
*Representative scales:* for LEO/VLEO, third-body perturbations are small compared with drag, but are systematically modelled for high-precision OD; JPL DE ephemerides describe alignment/accuracy relative to ICRF and data sources. citeturn5search1turn5search21  
*Correlation scales:* ephemeris errors are slowly varying (days–years) and highly correlated. citeturn5search21turn5search1  
*Entry into POD/propagation:* deterministic \(\mathbf{a}_b\) term with parametric uncertainty in \(\mathbf{r}_b(t)\). citeturn5search1

**Solar radiation pressure (SRP) and Earth albedo/IR uncertainty (including optical properties)**  
*Uncertainty model:* SRP acceleration in classic IERS formulation is
\[
\mathbf{a}_{SRP}=K\left(\frac{AU}{r_\odot}\right)^2 C_R \frac{A_\perp}{m}\,\hat{\mathbf{s}},
\quad K \approx 4.56\times 10^{-6}\,\mathrm{N/m^2}.
\]
\(C_R\) and effective area are uncertain and often estimated as a scale factor (constant per arc or RW/GM(1)). citeturn7search0turn7search14turn6search2  
Earth radiation pressure (ERP) uses Earth albedo+IR fluxes; modern LEO work increasingly uses CERES-driven flux products and simplified coefficient models, with uncertainties represented as scale factors on albedo/IR coefficients and/or time-varying stochastic residual accelerations. citeturn2search4turn2search7turn2search15  
*Representative magnitudes (VLEO):* SRP acceleration scales with \(A/m\) and is often \(10^{-8}\)–\(10^{-7}\,\mathrm{m/s^2}\) for typical LEO spacecraft, while ERP can be a significant fraction of SRP depending on geometry and altitude; in VLEO it is typically smaller than drag but relevant for high-precision residuals. citeturn7search0turn2search7turn2search15  
*Correlation scales:* SRP/ERP residuals show strong 1‑cycle‑per‑revolution signatures and eclipse-related structure; uncertainty is thus time-correlated on orbital timescales and often treated with empirical accelerations or RW coefficients. citeturn2search9turn12search3turn7search14  
*Entry into POD/propagation:* via \(\mathbf{a}_{RP}\) in dynamics and via torque modelling when attitude dynamics is estimated. citeturn2search15turn6search2

**Magnetic field model and residual dipole uncertainty**  
*Uncertainty model:* disturbance torque
\[
\mathbf{T}_{mag}=\mathbf{m}_{res}\times \mathbf{B}(\mathbf{r},t),
\]
with \(\mathbf{B}\) from a geomagnetic model (IGRF-like) and \(\mathbf{m}_{res}\) uncertain (pre-launch measurement error, in-orbit drift). \(\mathbf{m}_{res}\) is frequently modelled as constant plus RW/GM(1) drift. citeturn5search2turn5search10turn1search7  
*Representative ranges:* NASA’s classic magnetic torque report discusses disturbance torques from residual spacecraft dipoles; CubeSat-focused work similarly formulates torque scaling with residual dipole. citeturn5search2turn5search18 IGRF uncertainty estimation is addressed in recent peer-reviewed work quantifying uncertainty of the reference field coefficients. citeturn1search7  
*Correlation scales:* \(\mathbf{B}\) varies primarily with orbit position (periodic) plus storm-time external-field contributions; \(\mathbf{m}_{res}\) drift is slow (hours–months). citeturn1search7turn5search10  
*Entry into POD/propagation:* primarily attitude dynamics and pointing; indirectly affects POD when attitude uncertainty feeds into drag/RP modelling. citeturn5search2turn8search35

**Structural flexing and thermal deformation**  
*Uncertainty model:* typically expressed as time-varying misalignment between sensor frames and body frame, e.g., a small-angle vector \(\boldsymbol{\epsilon}(t)\) with RW/GM(1) drift:
\[
\mathbf{R}_{SB}(t)=\mathbf{R}_{SB}^0\exp([\boldsymbol{\epsilon}(t)]_\times),
\quad \boldsymbol{\epsilon}_{k+1}=\boldsymbol{\epsilon}_k+\eta_k.
\]
This is used for star-tracker alignment error, payload boresight drift, and appendage-induced attitude dynamics mismatch. citeturn4search11turn4search15turn4search7  
*Representative ranges:* ESA pointing error engineering examples treat alignment errors and their evolution as explicit error-budget terms, with drift-like behaviour addressed via calibration models. citeturn4search11turn4search15  
*Correlation scales:* thermally driven deformations correlate with orbital beta angle, eclipse transitions, and heater cycles (orbital period to seasonal). citeturn4search11turn4search15  
*Entry into POD/propagation:* through attitude determination errors (thus \(A_\perp\) in drag and incidence for RP), and through lever‑arm errors for GNSS antenna phase centre and accelerometer location (mapping into observation models). citeturn12search3turn4search11

**Mass property and fuel gauging uncertainty (mass \(m\), inertia \(\mathbf{I}\), propellant remaining)**  
*Uncertainty model:* mass enters as \(1/m\) in acceleration; inertia enters attitude dynamics. A standard treatment is:
\[
m(t)=m_{nom}(t)+\delta m(t),\quad \delta m\ \text{as bias or RW};
\qquad
\mathbf{I}(t)=\mathbf{I}_{nom}(t)+\delta\mathbf{I}.
\]
Fuel gauging uncertainty is often treated as a bounded bias (percent full-scale) updated intermittently, sometimes modelled as RW between updates. citeturn5search38turn13search11turn13search16  
*Representative ranges:* NASA propellant gauging technology summaries report uncertainties of order **~1–3%** in some contexts and note configuration dependence; modern microgravity gauging experiments report uncertainties like ±1% (1g) and larger (e.g., ±6% in low gravity) for particular methods. citeturn13search16turn13search11 Peer-reviewed gauging uncertainty analyses show Monte Carlo propagation of pressure/temperature sensor uncertainties into sub‑percent to percent-level residue volume uncertainty depending on method and instrumentation. citeturn13search0turn13search4  
*Correlation scales:* gauging bias is slow (days–months) with step changes at burns or slosh events; inertia uncertainty changes with propellant usage and slosh state. citeturn5search38turn5search34  
*Entry into POD/propagation:* directly in \(\mathbf{a}_{th}= \mathbf{F}/m\), \(\mathbf{a}_{drag}\propto 1/m\), and in attitude equations \(\mathbf{I}\dot{\boldsymbol{\omega}}+\boldsymbol{\omega}\times(\mathbf{I}\boldsymbol{\omega})=\mathbf{T}\). citeturn7search25turn8search35

**Propulsion thrust magnitude, direction, and thrust noise**  
*Uncertainty model:* finite burn modelling uses thrust \(\mathbf{F}(t)\) and mass flow \(\dot{m}\); thrust uncertainty is commonly decomposed into magnitude scale, pointing (misalignment), and noise:
\[
\mathbf{F}(t)=(1+\delta_F)\,\mathbf{R}(\boldsymbol{\delta\theta})\,\mathbf{F}_{cmd}(t)+\mathbf{n}_F(t),
\]
with \(\delta_F\) (dimensionless) as bias/RW, \(\boldsymbol{\delta\theta}\) (rad) as small-angle bias/RW, and \(\mathbf{n}_F\) as coloured noise if thrust ripple/noise matters. CCSDS navigation conventions provide the standard finite‑burn \(\Delta\mathbf{v}\) computation relationships used in OD interfaces. citeturn7search25turn7search1  
*Representative magnitudes:* the thrust domain depends strongly on propulsion type; peer-reviewed and standards-oriented literature on electric propulsion thrust measurement explicitly treats micro‑Newton to milli‑Newton thrust levels and emphasises the significance of thrust noise spectra when relevant. citeturn13search5turn13search25  
*Correlation scales:* thrust scale/misalignment drift is slow (thermal/mechanical); thrust noise can be broadband and sometimes characterised in frequency domain. citeturn13search5turn13search25  
*Entry into POD/propagation:* as \(\mathbf{a}_{th}=\mathbf{F}/m\) during burn, with event-time uncertainty handled as uncertain burn start/stop epochs; manoeuvre modelling choices (impulsive vs finite) matter for estimation. citeturn7search25turn4search34

**GNSS measurement errors (code/phase, multipath, receiver biases)**  
*Uncertainty model:* standard observation models (for satellite \(s\)) are
\[
P^s=\rho^s+c(\delta t_r-\delta t_s)+T^s+I^s+b_P+\epsilon_P,
\]
\[
\Phi^s=\rho^s+c(\delta t_r-\delta t_s)+T^s-I^s+\lambda N^s+b_\Phi+\epsilon_\Phi,
\]
with \(\epsilon\) often WGN, and biases \(b_P,b_\Phi\) including hardware delays and phase biases (often treated as constant or slow RW). citeturn4search24turn3search31turn3search11  
*Representative ranges:* GNSS measurement noise is typically described as **dm–m level for pseudorange** and **mm‑level for carrier phase** in LEO POD literature; multipath introduces systematic azimuth/elevation dependencies and can create non‑Gaussian tails. citeturn4search16turn4search28turn3search31 Example real-time POD work uses pseudorange and (epoch‑differenced) carrier-phase noise levels on the order of 1 m and 0.01 m respectively in algorithm descriptions, reflecting typical weighting assumptions. citeturn12search36turn4search36  
*Correlation scales:* multipath and uncalibrated antenna phase-centre variations are time-correlated with geometry (orbital period and attitude); ionospheric residuals after dual-frequency combinations can remain correlated spatially/temporally. citeturn4search24turn8search24turn3search7  
*Entry into POD/propagation:* through observation residuals and measurement covariance \(R_k\), often elevation-dependent. citeturn4search24turn4search36

**Satellite Laser Ranging (SLR) measurement errors**  
*Uncertainty model:* SLR range
\[
\rho_{SLR}=|\mathbf{r}-\mathbf{r}_{sta}|+\Delta_{tropo}+\Delta_{rel}+\Delta_{CoM}+b_{sta}+\epsilon,
\]
with \(\epsilon\) including single-shot timing noise and \(b_{sta}\) station-dependent biases; systematic detector/timing effects can create non-Gaussian residual structure. citeturn3search16turn3search20turn4search29  
*Representative ranges:* modern SLR systems discuss achieving **~mm normal-point precision** through averaging many shots; single-shot precision can be centimetres with mm-level normal points under high repetition and stable timing, but systematic errors remain limiting. citeturn3search16turn3search28  
*Correlation scales:* station biases can be stable over weeks–months but drift; detector artefacts and timing drifts can correlate over a pass. citeturn4search29turn3search20turn3search5  
*Entry into POD/propagation:* as high-precision range observations used for orbit validation and bias estimation; covariance \(R\) often separated into statistical + bias components. citeturn3search1turn4search21

**Star tracker, IMU, magnetometer sensor errors (bias, noise, scale factors)**  
*Uncertainty model:* common attitude sensor models:
- star tracker provides attitude or unit-vector measurements with WGN plus calibration bias;  
- gyros: \(\boldsymbol{\omega}_m=\boldsymbol{\omega}+ \mathbf{b}_g+\mathbf{n}_g\), with \(\dot{\mathbf{b}}_g=\mathbf{w}_b\) (RW) or GM(1);  
- magnetometer: \(\mathbf{B}_m=\mathbf{S}\mathbf{B}+\mathbf{b}_m+\mathbf{n}_m\), with scale/misalignment matrices \(\mathbf{S}\) and bias \(\mathbf{b}_m\). citeturn3search29turn3search17turn5search10  
Star tracker sampling rates of ~1–10 Hz and arcsecond-class attitude accuracy (in representative literature) motivate treating measurement noise as small-angle Gaussian, while misalignment drift is RW-like. citeturn3search29turn4search11  
*Representative parameter ranges:* IMU noise is often parameterised by angular random walk (ARW) and bias drift; Allan-variance-based methods are used to estimate conservative noise parameters, and RW bias is a standard assumption for many gyros. citeturn3search26turn3search38turn3search17  
Residual dipole and magnetometer calibration uncertainty are a recurring dominant term in magnetic-attitude error budgets for small spacecraft, motivating explicit bias/scale modelling. citeturn5search10turn5search3  
*Correlation scales:* gyro bias RW integrates over time (variance \(\propto t\)); magnetometer bias drift is slow; star tracker misalignment correlates with thermal state and can be effectively non-observable without manoeuvres. citeturn3search9turn4search15turn3search17  
*Entry into POD/propagation:* attitude uncertainty enters translational dynamics via attitude-dependent drag/RP area and in measurement models via antenna/sensor frame transformations. citeturn8search35turn4search11

**Timing and clock errors (receiver clock, time-tag errors)**  
*Uncertainty model:* clock bias \(b_c\) and drift \(\dot{b}_c\) are classically modelled as RW / integrated RW / Gauss‑Markov processes in navigation filtering:
\[
b_{c,k+1}=b_{c,k}+\dot{b}_{c,k}\Delta t + \eta_{b},\quad
\dot{b}_{c,k+1}=\dot{b}_{c,k}+\eta_{\dot{b}},
\]
with process noise defined per standard clock stochastic models. citeturn6search24turn3search11  
*Representative ranges:* IGS guidance documents discuss achieving high-precision clock estimates (order 0.1 ns or better in geodetic-grade contexts), which maps to centimetre-level range units via \(c\,\delta t\). citeturn3search11turn3search4  
*Correlation scales:* clock noise includes multiple power-law components; RW-type growth makes clock uncertainty increase with prediction horizon. citeturn6search24turn4search0  
*Entry into POD/propagation:* directly in GNSS \(P,\Phi\) via \(c(\delta t_r-\delta t_s)\), and indirectly through time-tag errors in SLR. citeturn4search24turn3search5

**Modelling errors and parameter estimation uncertainty (correlations, identifiability)**  
*Uncertainty model:* key OD literature identifies strong correlations between \(\rho\), \(C_D\), attitude, and empirical accelerations. One formulation is joint estimation of density and drag-coefficient for time-varying attitude, using augmented-state filters and careful modelling of parameter processes (RW/GM(1)). citeturn12search19turn8search35turn6search19  
For covariance realism, recent work introduces stochastic consider-parameter techniques and explicit time-correlated error modelling to avoid underestimating prediction uncertainty when parameters are treated as deterministic. citeturn12search5turn12search16turn12search32  
*Representative quantitative statements:* scale differences between thermospheric datasets are attributed in part to aerodynamic modelling (geometry + GSI), implying that “effective” density errors can reflect \(C_D\) errors; this is treated via stochastic parameterisation rather than pure measurement noise. citeturn8search35turn12search17  
*Correlation scales:* identifiability depends on tracking cadence and attitude excitation; time-correlated process modelling (with uncertain \(\tau\)) is explicitly studied as a driver of (mis)estimated confidence bounds in Kalman filtering. citeturn12search32turn6search13  
*Entry into POD/propagation:* appears as augmented parameters in \(\boldsymbol{\theta}\) and in the structure of \(\mathbf{Q}\) and \(\mathbf{R}\) used in the estimator. citeturn6search13turn12search5

**Operational uncertainties (manoeuvre execution, ground coverage, tracking geometry, latency)**  
*Uncertainty model:* operational uncertainties are often represented as (i) impulsive/finite-burn \(\Delta\mathbf{v}\) execution error with covariance in magnitude/direction/time, (ii) measurement scheduling as a random availability process (Bernoulli dropout), and (iii) latency as delayed measurements in filtering.
CCSDS Orbit Data Messages explicitly support covariance exchange and manoeuvre specification; this formalises how operationally derived uncertainty is packaged and propagated. citeturn7search1turn7search25  
Tracking geometry effects are treated via observability/conditioning: measurement noise maps to state covariance depending on pass geometry, measurement types, and data gaps; OD error analyses study how station passes and timeliness constraints affect orbit refinement. citeturn12search28turn7search22  
*Representative ranges:* manoeuvre execution uncertainty is mission/actuator dependent; the literature treats it through Δv covariance and through detection/estimation methods for unknown manoeuvres. citeturn4search18turn7search25  
*Correlation scales:* outage processes can be bursty (storm-time GNSS, comms); latency introduces effective time correlation in residuals when not properly handled. citeturn12search36turn7search16  
*Entry into POD/propagation:* as uncertain state discontinuities \(\Delta\mathbf{v}\) and as time-varying measurement operator sets \(\mathbf{H}_k\) and covariances \(R_k\). citeturn7search25turn12search28

**Software and human factors (statistical descriptions where available)**  
*Uncertainty model:* software/human contributions are often modelled as **rare-event risk** and categorical failure modes rather than Gaussian noise; statistically, they are captured by (i) failure-rate models, (ii) discrete-event processes causing state/configuration jumps, and (iii) heavy-tailed distributions for error magnitude when incidents occur. citeturn7search21turn7search3turn7search6  
*Representative quantitative statements:* reviews of spaceflight failures report that software is responsible for a non-trivial fraction of failures (with reported ranges across studies), while NASA human-error analysis work emphasises systematic contributors and mishap statistics in aerospace contexts. citeturn7search18turn7search15  
*Correlation scales:* such errors are not stationary; they are event-driven and clustered around operationally stressful intervals (critical manoeuvres, configuration changes, anomalies). citeturn7search3turn7search21  
*Entry into POD/propagation:* through wrong/late commands (wrong thrust, wrong attitude mode), wrong ground processing parameters (bias injection), and data-handling faults (time-tag errors, frame mismatches), typically represented as discrete changes in \(\boldsymbol{\theta}\) or state. citeturn7search3turn7search21

**Measurement biases and stochastic noise (generic cross-cutting entry)**  
*Uncertainty model:* many systems adopt a unified measurement error model
\[
\mathbf{y}=\mathbf{h}(\mathbf{x})+\mathbf{b}+\boldsymbol{\nu},\quad
\mathbf{b}_{k+1}=\mathbf{b}_k+\eta_k,
\]
with \(\boldsymbol{\nu}\) Gaussian and \(\eta_k\) RW; for robustness, \(\boldsymbol{\nu}\) is sometimes modelled as Student‑t or as a Gaussian mixture to capture outliers. Time-correlated measurement/process errors are explicitly addressed in modern estimation bounds work for Gauss–Markov structures with uncertain correlation times. citeturn12search32turn7search29turn3search20  
*Representative ranges:* bias magnitudes depend on sensor; SLR range biases at stations are assessed at mm–cm scales, GNSS biases are represented via observable-specific bias products (OSB) and receiver hardware delays, and unmodelled accelerations are handled via stochastic processes. citeturn4search29turn3search31turn6search13  
*Correlation scales:* biases can persist for weeks–months; noise can be white at high frequency but becomes coloured after differencing, filtering, and due to geometry-driven systematic effects. citeturn12search3turn4search16turn3search20  
*Entry into POD/propagation:* through \(R_k\), through explicit bias states, and through empirical accelerations / process noise shaping in \(Q_k\). citeturn7search29turn6search13

## Worked quantitative examples

### Density scale factor as multiplicative bias and its effect on drag acceleration

Assume a simple drag model with \(A_\perp/m = 0.02\,\mathrm{m^2/kg}\), \(C_D=2.2\), \(v_{rel}=7.7\,\mathrm{km/s}\). The drag acceleration magnitude is
\[
a_D=\tfrac12 \rho\,C_D\,(A/m)\,v_{rel}^2.
\]
At ~400–450 km, the NASA profile figure implies densities on the order of \(10^{-13}\)–\(10^{-12}\,\mathrm{kg/m^3}\) depending on solar activity and diurnal state. citeturn11view0  
Take \(\rho=5\times10^{-13}\,\mathrm{kg/m^3}\). Then
\[
a_D \approx \tfrac12 (5\!\times\!10^{-13})\cdot 2.2\cdot 0.02\cdot (7.7\!\times\!10^3)^2
\approx 6.5\times10^{-7}\,\mathrm{m/s^2}.
\]
If density is modelled with a multiplicative scale \(s_\rho=\exp(\delta_\rho)\) and \(\delta_\rho\) has 1‑σ of 0.2 (≈20% in small-error approximation), then \(a_D\) inherits ≈20% 1‑σ uncertainty, i.e. \(\sigma_{a_D}\approx 1.3\times10^{-7}\,\mathrm{m/s^2}\). Such multiplicative density uncertainty levels (quiet ~10–15% or larger) are routinely cited for thermosphere modelling. citeturn12search18turn12search26

### Drag coefficient \(C_D\) as a random walk parameter in OD

A common fitted-parameter OD model is
\[
C_{D,k+1}=C_{D,k}+\eta_k,\quad \eta_k\sim \mathcal{N}(0,q_{C_D}\Delta t).
\]
Suppose \(q_{C_D}=10^{-6}\,\mathrm{s^{-1}}\) (illustrative tuned value) and \(\Delta t=600\,\mathrm{s}\). Then \(\sigma_{\Delta C_D}=\sqrt{q\Delta t}\approx 0.024\) per 10 minutes. Over \(N\) steps, RW variance grows linearly: \(\mathrm{Var}(C_D)\approx N q\Delta t\). This RW modelling strategy is explicitly documented and contrasted in real-time POD configurations where drag coefficient is treated as RW. citeturn7search14turn7search5turn6search19

### GNSS carrier-phase noise mapped to range/position uncertainty

Carrier phase observation noise \(\sigma_\Phi\) propagates roughly to line-of-sight range noise \(\sigma_\rho\approx \sigma_\Phi\) (in metres) for the phase observable (before ambiguity handling). Many LEO real-time OD algorithms weight carrier phase much higher than code, e.g. adopting \(\sigma_P\sim 1\,\mathrm{m}\) and \(\sigma_{\Delta\Phi}\sim 0.01\,\mathrm{m}\) for epoch-differenced phase in simplified real-time filters. citeturn12search36turn4search36  
If geometry is good, a back-of-envelope position error scale is
\[
\sigma_{pos}\approx \frac{\sigma_\rho}{\sqrt{N_{sat}}}\cdot \mathrm{GDOP},
\]
so with \(N_{sat}=8\), \(\sigma_\rho=0.01\,\mathrm{m}\), and GDOP ≈ 2, \(\sigma_{pos}\sim 7\,\mathrm{mm}\) (measurement-only). In practice for VLEO, dynamic model errors (drag) dominate unless statistically modelled/estimated appropriately. citeturn4search28turn7search29turn12search18

### Clock bias error to range error and why clock processes are modelled as RW/GM

In GNSS, a receiver clock bias \(\delta t\) contributes \(c\,\delta t\) to pseudorange. A 1 ns bias corresponds to ~0.30 m range. Clock models used in navigation commonly represent clock bias/drift as RW/integrated RW or Gauss–Markov processes; explicit forms and their covariance growth are documented in standard treatments of clock stochastic processes. citeturn6search24turn3search11  
This is why many OD filters include clock states with process noise, rather than treating clock offset as an independent WGN term. citeturn3search11turn6search24

## Comparative summary table

The table uses compressed phraseology; ranges are “representative” of values reported or implied in the cited sources and should be interpreted as order-of-magnitude guidance for 150–450 km.

| Source | Typical stochastic model type | Typical parameter ranges (units) | Dominant timescale | Spatial correlation scale | Key references |
|---|---|---|---|---|---|
| Neutral density \(\rho\) | multiplicative scale \(s_\rho=\exp(\delta)\), \(\delta\) as RW/GM(1) | \(\rho\sim10^{-9}\) kg m⁻³ @ ~200 km; \(\rho\sim10^{-13}\)–\(10^{-12}\) kg m⁻³ @ ~400–450 km; quiet-model error ~10–15% | minutes–days; diurnal & storm | anisotropic; along-track aggregation ~10⁴ km in some retrieval uncertainty analyses | citeturn11view0turn12search18turn12search3 |
| Winds \(\mathbf{u}\) | additive WGN/GM(1) on wind components; often absorbed into drag residual | “hundreds m/s” max velocities noted; model–obs differences can be large seasonally | minutes–hours (waves), hours (tides) | regional; global patterns in empirical wind models | citeturn0search1turn9view0turn0search24 |
| Composition (\(\Sigma O/N_2\), \(\bar m\)) | correlated deviations from climatology; often folded into \(s_\rho\) and \(C_D\) uncertainty | \(\Sigma O/N_2\) order ~0.4–0.8 (context-dependent) | hours–days (disturbance transport) | global / hemispheric during storms | citeturn5search4turn5search0turn5search35 |
| \(C_D\), accommodation, scattering kernel | \(C_D\) RW/GM(1); kernel-model selection; uncertain accommodation params | \(C_D\sim1.2\)–2.5 (geometry/altitude dependent); strong sensitivity to specularity | hours–months (surface state), seconds (attitude) | local to spacecraft; correlated with attitude | citeturn8search35turn7search5turn8search2 |
| Rarefied/Kn regime | model-form uncertainty; uncertain \(Kn\) via uncertain \(\lambda\) | \(Kn\) explicitly used to classify VLEO flow regime | altitude-driven; modulated by density | altitude/local-time dependent | citeturn8search25turn8search33turn8search2 |
| Atomic oxygen flux/erosion | flux via \(n\,v\); erosion via yield×fluence; multiplicative uncertainty | AO flux ~10¹⁴–10¹⁵ atoms cm⁻² s⁻¹ (LEO) | minutes–years (fluence accumulation) | follows thermosphere; global patterns | citeturn2search6turn2search3turn2search25 |
| Charging/plasma | regime/event-driven; bounded uncertain plasma params; event processes | charging concerns noted for energetic auroral particles; VLEO usually lower potentials than high altitudes | seconds–minutes; eclipse/storm | auroral regions; local-time dependence | citeturn9view2turn0search25turn9view0 |
| Aerodynamic residual accel/torque | empirical accel as GM(1), direction-dependent \(\tau\) | correlation times e.g. ~5400 s along-track, ~360 s cross/radial in one study | minutes–orbit | corresponds to orbital-path scales | citeturn8search20turn6search13turn7search29 |
| F10.7, Kp, Ap drivers | additive forecast error; indices piecewise constant (Kp/ap) | Kp 3‑hourly; Ap daily in 2 nT units; F10.7 in sfu | 3 h bins (Kp/ap), days (F10.7) | global driver values | citeturn1search1turn1search22turn1search9 |
| Storm-time response | jump / regime-switch bias on log-density; heavy tails | global density response can be O(10–100%) or larger depending on event | hours–days | global/hemispheric | citeturn9view0turn6search12turn6search7 |
| Geopotential higher orders | coefficient covariance; linearised \(\delta C_{\ell m}\) | EGM2008 provides high-degree model with error information | static (model), periodic orbit signatures | wavelength \(\sim 2\pi R/\ell\) | citeturn1search0turn1search28turn1search12 |
| Earth tides/loading | deterministic periodic + residual coloured error | cm–dm displacements; mm-level residual relevance for SLR/GNSS | tidal periods; days (loading) | global, site-dependent | citeturn0search26turn0search14turn3search1 |
| Third-body ephemerides | parametric uncertainty in \(\mathbf{r}_b(t)\) | DE ephemerides provide accuracy context | days–years | global (solar system) | citeturn5search1turn5search21 |
| SRP | scale factor \(C_R\) RW/GM(1); coloured residual | \(K=4.56\times10^{-6}\,\mathrm{N/m^2}\); accel scales with \(A/m\) | orbit period; eclipse | global | citeturn7search0turn7search14turn6search2 |
| Earth albedo/IR | flux-model uncertainty + coefficient scale factors | CERES-based modelling; simplified coefficient uncertainty | orbital+seasonal | hemispheric/regional | citeturn2search4turn2search7turn2search15 |
| Magnetic field & residual dipole | \(\mathbf{m}_{res}\) bias/RW; \(\mathbf{B}\) model uncertainty | IGRF coefficient uncertainty; residual dipole dominates small-sat torques | orbit-period + months | global; storm-time external fields | citeturn1search7turn5search2turn5search10 |
| Structural/thermal deformation | misalignment RW/GM(1) | alignment error terms treated explicitly in pointing budgets | orbit-period to seasonal | spacecraft-specific | citeturn4search11turn4search15turn4search7 |
| Mass & fuel gauging | bias + step/RW at burns | ~1–3% (context-dependent); method-specific values | days–months | N/A | citeturn13search16turn13search11turn13search0 |
| Thrust magnitude/direction/noise | scale + misalignment + coloured noise | micro‑N to mN thrust regimes for EP; thrust noise spectrum relevant | seconds–hours (noise), long-term drift | N/A | citeturn7search25turn13search5turn13search25 |
| GNSS code/phase/multipath/bias | WGN + geometry-correlated systematics; biases constant/RW | code dm–m; phase mm‑level; strong multipath systematic patterns | seconds–orbit | geometry-driven | citeturn4search24turn4search16turn3search31 |
| SLR errors | shot noise + station bias + detector artefacts | mm normal points possible; station RBs mm–cm | pass–months | station-specific | citeturn3search16turn3search20turn4search29 |
| Star tracker | small-angle WGN + misalignment RW | arcsecond-class accuracy typical; 1–10 Hz sampling | seconds–orbit; thermal | spacecraft-specific | citeturn3search29turn4search11 |
| IMU (gyro/accel) | WGN + bias RW/GM(1); Allan-variance parameterisation | ARW & bias drift parameters; RW bias common | seconds–hours | N/A | citeturn3search26turn3search38turn3search17 |
| Magnetometer | scale/bias/misalignment; residual dipole coupling | calibration/magnetic cleanliness dominates | orbit–months | local geomagnetic | citeturn5search10turn5search2turn1search7 |
| Clock/time-tag | RW/IRW/GM; bias+drift states | 1 ns → 0.30 m range; high-grade clocks much better | seconds–days | N/A | citeturn6search24turn3search11turn3search5 |
| Estimation identifiability | parameter correlations; stochastic consider parameters | strong \(\rho\)–\(C_D\) correlation; time-correlation uncertainty | arc length dependent | N/A | citeturn12search19turn12search5turn12search32 |
| Ops (coverage/latency) | dropout processes; Δv timing/mag covariance | geometry affects covariance; CDM-based statistics used | minutes–days | network-dependent | citeturn7search1turn12search28turn7search16 |
| Software & human | rare-event / categorical; heavy-tailed consequences | software 3–33% of failures (range across studies); human-error mishap emphasis | event-driven | N/A | citeturn7search18turn7search3turn7search15 |

## Notes on evidence quality and practical interpretation

The most consistently quantified uncertainty components in public literature for LEO/VLEO OD are: (i) **thermosphere model error statistics** under quiet/storm regimes, (ii) **stochastic OD process-noise architectures** (RW/GM(1)) for empirical accelerations and coefficients, and (iii) **tracking measurement noise/bias models** for GNSS/SLR. citeturn12search18turn6search13turn3search20

By contrast, **spatial correlation lengths** for neutral-density *errors* (not just variability) and for GSI parameter drift are less standardised: some studies report effective along-track resolutions/averaging scales (e.g., ~12,000 km) or direction-dependent correlation times (e.g., one orbit along-track), but these often blend physics with retrieval/processing noise. citeturn12search3turn8search20turn12search29

Finally, several uncertainties of high operational importance—particularly **software/human factors** and **rare failure modes**—do not naturally fit a Gaussian small-noise paradigm. The literature typically frames them as categorical incident classes and reliability statistics, motivating discrete-event and heavy-tailed modelling if one wants to integrate them into end-to-end probabilistic OD/control risk models. citeturn7search18turn7search3turn7search21