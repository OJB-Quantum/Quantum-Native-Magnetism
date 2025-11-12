# Quantum-Native-Magnetism
A repository for performing quantum-native scientific computing for atomistic and nanoscale magnetism studies. Created by Onri Jay Benally. 

---

The following techniques are at least eligible for quantum-native magnetic modeling on Noisy Intermediate-Scale Quantum (NISQ) computing systems:

```
- Thermofield-Double State (purified Gibbs state, VQE and VQS-based, NISQ and Near-Term Quantum) 
- Quantum Imaginary-Time Evolution (QITE) (finite-temperature/ground-state prep, NISQ) 
- Stochastic Schrödinger‑equation discretizations for a few modes (QSDE prototype, NISQ) 
- Quantum Trajectory (NISQ and Near-Term Quantum) 
- Short-Depth Trotter (NISQ and Near-Term Quantum) 
```

One of the best techniques to implement for NISQ systems can be achieved with the following:
- QITE implemented with error mitigation on as many qubits as possible.


---

The following techniques are eligible for quantum-native magnetic modeling on Near-Term or Fault-Tolerant quantum computing systems: 

```
- Stochastic Dirac‑equation sampling (Near-Term Quantum) 
- Quantum Metropolis Sampling (preparing thermal distributions, for Near-Term) 
- Quantum stochastic differential equations (QSDE) (Fault-Tolerant Quantum)
```

---



| Aspect | Ehrenfest–LLB–Boltzmann (E‑LLB‑B) | Ehrenfest–LL–Boltzmann (E‑LL‑B) |
|---|---|---|
| Drift structure | Precession + transverse relaxation ($\Gamma_\perp$) + **longitudinal relaxation ($\Gamma_\parallel$)** toward $m_eq(T)$ | Precession + transverse relaxation ($\Gamma_\perp$) only; **no longitudinal channel** |
| Stationary thermodynamics | Gibbs/Boltzmann fixed point when rates obey detailed balance (KMS) | Orientation-only equilibrium on a fixed-magnitude shell; no intrinsic route to m -> 0 |
| Critical regime near $T_c$ | Captures collapse of magnetization magnitude; supports $\chi(T)$ peaks and Binder cumulants | Misses magnitude collapse; cannot reproduce Tc without adding a longitudinal term |
| Kinetic inputs | Can include magnon/electron/phonon collision integrals; material-specific rates | Can include collisions, but without a longitudinal channel critical behavior is still absent |
| Use case | Curie-temperature prediction, finite-size scaling, non-equilibrium transients | Directional diffusion and transport away from criticality; **not** reliable for $T_c$ |

---

| Model | Longitudinal relaxation present? | Guaranteed thermal fixed point (Gibbs/Boltzmann)? | Finite-size critical diagnostics ($\chi$, Binder) | Non-equilibrium kinetics (carriers/magnons) | $T_c$ capability | Overall comprehensiveness |
|---|---|---|---|---|---|---|
| **Ehrenfest–LL–Boltzmann (E‑LL‑B)** | No (orientation-only) | Not generically near $T_c$ (no $m_{eq}(T)$ channel) | Partial (orientation diffusion only) | Limited | **Weak near $T_c$** unless upgraded to include a longitudinal channel | Narrow; misses magnitude collapse that defines Tc |
| **Kinetic‑aware Ehrenfest–LLB–Boltzmann (E‑LLB‑B)** | **Yes** ($\Gamma_\parallel(T)$) | **Yes** with detailed-balance Lindbladian anchoring | **Yes** (full $P(m)$ dynamics enables $\chi$ and Binder) | **Yes** via collision-integral-based rates (qLLB-style) | **Strong** (robust $T_c$ from $m(T)$, $\chi(T)$, Binder $U_4$) | Broad: equilibrium + kinetics + finite-size scaling |
| **Ehrenfest–LLB (GKSL, expectation-level)** | **Yes** (via Lindblad rates) | **Yes** (KMS/Davies choice) | Moments only; sample trajectories for chi and Binder | Indirect (time-dependent rates) | **Good** if rates are calibrated from microscopic models | Medium–high: efficient but does not evolve full $P(m)$ |

---

Relevant equation: 

## Notation

* Magnetization (unit): ( $\boldsymbol{m}=\boldsymbol{M}/M_s$ ), ( $m=\lvert\boldsymbol{m}\rvert$ ).
* Effective field: ( $\boldsymbol{B}_{\mathrm{eff}}=-\partial \mathcal{F}/\partial \boldsymbol{m}$ ), with Zeeman, exchange, anisotropy, Dzyaloshinskii–Moriya (DMI), etc.
* Quantum expectation (site (i)): ( $\boldsymbol{s}_i(t)$ = $\mathrm{Tr}!\big[\rho(t),\boldsymbol{\sigma}_i\big]$ ).
* Unit field direction: ( $\hat{\boldsymbol{b}}$ = $\boldsymbol{B}*{\mathrm{eff}}$ / $\lvert \boldsymbol{B} * {\mathrm{eff}} \rvert$ ).
* Rates: longitudinal ( $\Gamma_\parallel=1/T_1$ ), transverse ( $\Gamma_\perp=1/T_2$ ); LLB coefficients ( $\alpha_\parallel(T),\alpha_\perp(T)$ ).
* Equilibrium magnetization: ( $m_{\mathrm{eq}}(T)$ ) (or ( $m_{\mathrm{eq}}(T,D)$ ) if a finite‑size parameter ($D$) is included).

---

## 1) Landau–Lifshitz/ Gilbert and stochastic LLG

### 1.1 Classical original (Gilbert and LL forms)

$$
\frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}} +\frac{\alpha}{m_s},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt} \quad \text{(Gilbert)},,
$$

$$
\frac{d\boldsymbol{m}}{dt} = -\gamma',\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} -\gamma'\alpha,\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}}\big), 
\quad \gamma'=\frac{\gamma}{1+\alpha^2}\quad\text{(LL)},. 
$$

### 1.2 Stochastic LLG (thermal field)

$$
\frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\big(\boldsymbol{B}*{\mathrm{eff}}+\boldsymbol{B}*{\mathrm{th}}(t)\big) +\frac{\alpha}{m_s},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt},
$$

$$
\big\langle B_{\mathrm{th},\mu}(t),B_{\mathrm{th},\nu}(t')\big\rangle =2D,\delta_{\mu\nu},\delta(t-t'),,\qquad D \propto \frac{\alpha k_{\mathrm B}T}{\gamma M_s V},.
$$

### 1.3 Ehrenfest (expectation) form for LL‑type transverse damping

$$
\frac{d\boldsymbol{s}}{dt} =\frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} -\Gamma*{\perp}\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})\hat{\boldsymbol{b}}\Big),.
$$

### 1.4 Adapted (Qiskit‑ready)

**Hamiltonian (single macrospin or per site (i))**

$$
H(t)= -\frac{\hbar\gamma}{2},\boldsymbol{B}_{\mathrm{eff}}(t)\cdot\boldsymbol{\sigma} \quad\text{(encode as a sum of Pauli terms)}.
$$

**Optional LL‑like transverse damping on simulators**
Use dephasing Lindblad $(L^z=\sqrt{\gamma_\phi},\sigma^z)$ to reproduce a chosen ($T_2$), see §4.

---

## 2) Landau–Lifshitz–Bloch (LLB)

### 2.1 Classical original (single‑sublattice macrospin)

$$
\frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}} * \frac{\gamma,\alpha_{\parallel}(T)}{m^2},(\boldsymbol{m}!\cdot!\boldsymbol{B}_{\mathrm{eff}}),\boldsymbol{m} - \frac{\gamma,\alpha_{\perp}(T)}{m^2},\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}\big).
$$

A standard longitudinal field used inside ( $\boldsymbol{B}*{\mathrm{eff}}$ ) is

$$
\boldsymbol{B}*{\parallel} =\frac{1}{\chi_{\parallel}(T)} \left(1-\frac{m^2}{m_{\mathrm{eq}}^2(T)}\right)\boldsymbol{m}, \quad\text{so}\quad \lvert \boldsymbol{m}\rvert \to m_{\mathrm{eq}}(T).
$$

### 2.2 **Ehrenfest form of LLB** (deterministic mean)

$$
\frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} -\Gamma*{\perp}(T)\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})\hat{\boldsymbol{b}}\Big) -\Gamma_{\parallel}(T)\Big[(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}},.
$$

### 2.3 Adapted (Qiskit‑ready)

**Jump operators and rates**

$$
L^{-}=\sqrt{\gamma_{\downarrow}},\sigma^{-},\qquad L^{+}=\sqrt{\gamma_{\uparrow}},\sigma^{+},\qquad L^{z}=\sqrt{\gamma_{\phi}},\sigma^{z},
$$

$$
\Gamma_{\parallel}=\gamma_{\downarrow}+\gamma_{\uparrow},\qquad \Gamma_{\perp}=\frac{\Gamma_{\parallel}}{2}+\gamma_{\phi},\qquad \frac{\gamma_{\uparrow}}{\gamma_{\downarrow}}=e^{-\beta\hbar\omega}\ \ (\text{KMS}).
$$

Implement with `qiskit_dynamics.models.LindbladModel` and a Pauli‑decomposed ($H$).

### 2.4 **Quantum‑native LLB (qLLB)** — deterministic mean with quantum‑derived rates

**Same drift as 2.2**, but now use **qLLB coefficients**; e.g., for a spin–phonon case:

$$
\alpha_{\parallel}(T) =\lambda,\frac{2T}{3T_C},\frac{2q_s}{\sinh(2q_s)},\qquad \alpha_{\perp}(T) =\lambda\left[\frac{\tanh q_s}{q_s}-\frac{2T}{3T_C}\right],\qquad q_s=\frac{\mu H_{\mathrm{MFA}}}{k_{\mathrm B}T},,
$$

with ($m_{\mathrm{eq}}(T)$) from the Brillouin/Langevin mean‑field appropriate to the material.

**Stochastic extension (qLLB‑SDE, trajectory level)**
Write anisotropic noise consistent with fluctuation–dissipation:

$$
d\boldsymbol{s} =\Big[\text{qLLB drift of 2.2 with qLLB rates}\Big],dt +\sqrt{2D_{\perp}(T)},(\boldsymbol{I}-\hat{\boldsymbol{b}}\hat{\boldsymbol{b}}^{!\top}),d\boldsymbol{W}*{\perp} +\sqrt{2D*{\parallel}(T)} \hat{\boldsymbol{b}},dW_{\parallel},
$$

with ( $\langle dW_\mu\rangle=0$ ), ( $\langle dW_\mu(t),dW_\nu(t)\rangle=\delta_{\mu\nu},dt$ ), and
( $D_{\parallel,\perp}(T)$ ) chosen so that the stationary distribution is Boltzmann for the LLB free energy (i.e., Fokker–Planck consistency).

---

## 3) Kinetic‑aware Ehrenfest models on the Bloch ball ( $f$($\boldsymbol{m}$, $t$) )

Define a drift–diffusion–collision kinetic equation

$$
\partial_t f +\nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{A}(\boldsymbol{m},T),f\big] =\nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{D}(\boldsymbol{m},T),\nabla_{\boldsymbol{m}} f\big] +\mathcal{C}[f],
$$

where ($\mathcal{C}[f]$) encodes magnon/electron/phonon collision integrals, and ( $\boldsymbol{D}$ ) satisfies fluctuation–dissipation.

### 3.1 **Ehrenfest–LLB–Boltzmann (E‑LLB‑B)**

Use **LLB drift**:

$$
\boldsymbol{A}*{\mathrm{LLB}}(\boldsymbol{m},T) =\frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} -\Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot! \hat{\boldsymbol{b}})\hat{\boldsymbol{b}}\Big) -\Gamma_{\parallel}(T)\Big[(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}},
$$

and let ($\mathcal{C}[f]$) (e.g., Orbach/Raman/direct spin‑lattice processes; magnon–phonon kernels) determine temperature‑ and material‑dependent rates (qLLB‑style).

### 3.2 **Ehrenfest–LL–Boltzmann (E‑LL‑B)**

Use **LL transverse drift only**:

$$
\boldsymbol{A}*{\mathrm{LL}}(\boldsymbol{m},T) =\frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} -\Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot! \hat{\boldsymbol{b}})\hat{\boldsymbol{b}}\Big)!,
$$

which lacks a longitudinal channel, so it cannot by itself produce the amplitude collapse ( $m\to 0$ ) at ( $T_C$ ).

---

## 4) Lindblad (GKSL) backbone and thermal anchoring

### 4.1 Quantum‑native GKSL master equation

$$
\dot{\rho} = -\frac{i}{\hbar}[H,\rho] * \sum_\mu \Big(L_\mu,\rho,L_\mu^\dagger-\tfrac12{L_\mu^\dagger L_\mu,\rho}\Big) \equiv -\frac{i}{\hbar}[H,\rho]+\sum_\mu \mathcal{D}[L_\mu]\rho.
$$

### 4.2 Detailed balance (KMS) and the ($T_1/T_2$) dictionary

For Bohr frequency ($\omega$) between system levels,

$$
\frac{\gamma_{\uparrow}}{\gamma_{\downarrow}}=e^{-\beta\hbar\omega},\quad \Gamma_{\parallel}=\gamma_{\downarrow}+\gamma_{\uparrow},\quad \Gamma_{\perp}=\frac{\Gamma_{\parallel}}{2}+\gamma_{\phi}.
$$

### 4.3 “LL from Lindblad,” and a quantum LLG alternative

* Properly chosen dissipators yield **LL‑type transverse drift** at the expectation level (bridge to §1.3).
* A **quantum LLG** master equation (purity‑preserving) exists as an alternative to GKSL when one wants LL‑type damping with conserved purity.

---

## 5) Adapted Hamiltonian building blocks (for Qiskit encodings)

Nearest‑neighbor Heisenberg + Zeeman + DMI on a qubit graph:

$$
H = -\frac{\hbar\gamma}{2}\sum_i \boldsymbol{B}*{\mathrm{eff},i}!\cdot!\boldsymbol{\sigma}*i -\sum*{\langle i,j\rangle}!\big(J_x X_iX_j+J_y Y_iY_j+J_z Z_iZ_j\big) +\sum*{\langle i,j\rangle}!\boldsymbol{D}_{ij}!\cdot!(\boldsymbol{\sigma}_i\times\boldsymbol{\sigma}_j) ; (+,\text{penalties/anisotropy}).
$$

**Lindblad channels (per site ($i$))**

$$
L_i^{-}=\sqrt{\gamma_{\downarrow,i}},\sigma_i^{-},\quad L_i^{+}=\sqrt{\gamma_{\uparrow,i}},\sigma_i^{+},\quad L_i^{z}=\sqrt{\gamma_{\phi,i}},\sigma_i^{z}.
$$

These three primitives (Pauli‑decomposed ($H$); ($L_i^\pm,L_i^z$) with KMS rates) reproduce the **Ehrenfest–LLB** drift in §2.2. Adding a kinetic outer loop for ($\mathcal{C}[f]$) turns it into **E‑LLB‑B** (§3.1).

---

## 6) Index of equations

* **LLG/LL (classical):** two forms in §1.1.
* **sLLG:** thermal‑field form + white‑noise correlator in §1.2.
* **Ehrenfest LL (transverse‑only):** §1.3.
* **LLB (classical):** §2.1 (with longitudinal field piece).
* **Ehrenfest LLB (deterministic mean):** §2.2.
* **Ehrenfest–qLLB (deterministic mean):** same as §2.2 with **qLLB** coefficients in §2.4.
* **Ehrenfest–qLLB (stochastic extension):** SDE in §2.4.
* **E‑LLB‑B kinetic PDE:** §3.1.
* **E‑LL‑B kinetic PDE (LL drift only):** §3.2.
* **GKSL (quantum‑native) + KMS:** §4.1–§4.2.
* **Adapted Hamiltonian + dissipators (Qiskit):** §5. 

---

### Symbol mini‑glossary

( $\gamma$ ): gyromagnetic ratio; ( $\alpha$ ): Gilbert damping; ( $\alpha_{\parallel,\perp}(T)$ ): LLB coefficients; ( $\Gamma_{\parallel,\perp}$ ): Ehrenfest rates; ( $m_{\mathrm{eq}}(T)$ ): equilibrium magnetization; ( $\chi_{\parallel}(T)$ ): longitudinal susceptibility; ( $\boldsymbol{\sigma}$ = $(X,Y,Z)$ ): Pauli vector; ( $J_{x,y,z}$ ): exchange couplings; ( $\boldsymbol{D} * {ij}$ ): DMI vector; ( $\gamma*{\uparrow,\downarrow,\phi}$ ): Lindblad rate parameters; ( $\beta=(k_{\mathrm B}T)^{-1}$ ).
