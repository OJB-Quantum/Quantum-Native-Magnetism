# Quantum-Native-Magnetism
A repository for performing scientific quantum computing for atomistic and nanoscale magnetism studies. Created by Onri Jay Benally. 

[![License: GNU GPL v3](https://img.shields.io/badge/GNU_General_Public-License-Green)](https://choosealicense.com/licenses/gpl-3.0/)


---

The following techniques are eligible for quantum-native magnetic modeling on Noisy Intermediate-Scale Quantum (NISQ) computing systems:

```
- Thermofield-Double State (purified Gibbs state, VQE and VQS-based, NISQ and Near-Term Quantum) 
- Quantum Imaginary-Time Evolution (QITE) (finite-temperature/ground-state prep, NISQ) 
- Stochastic Schrödinger‑equation discretizations for a few modes (QSDE prototype, NISQ) 
- Quantum Trajectory (NISQ and Near-Term Quantum) 
- Short-Depth Trotter (NISQ and Near-Term Quantum) 
```

Four of the best techniques to implement for NISQ systems can be achieved with the following: 

```
- Quantum Imaginary-Time Evolution (QITE) on a NISQ quantum computer implemented with quantum error mitigation (QEM) on as many qubits as possible.
- Thermo Field Double (TFD) State on a NISQ or near-term quantum computer implemented with quantum error mitigation (QEM) or quantum error correction (QEC).
- Stochastic Schrödinger‑Equation-(SSE)-Discretization-implemented near-term quantum computer with quantum error mitigation (QEM) or quantum error correction (QEC).
- Quantum-Metropolis-Sampling-(QMS)-implemented near-term quantum computer with quantum error mitigation (QEM) or quantum error correction (QEC).
```

The following techniques are eligible for quantum-native magnetic modeling on Near-Term or Fault-Tolerant quantum computing systems: 

```
- Stochastic Dirac‑equation sampling (Near-Term Quantum) 
- Quantum Metropolis Sampling (preparing thermal distributions, for Near-Term) 
- Quantum stochastic differential equations (QSDE) (Fault-Tolerant Quantum)
```


* **Near-term:**
  Ehrenfest–LLB–Boltzmann + QITE + QEM on a **156-qubit Heron-class processor** (for example, `ibm_fez` or its siblings in the Heron r2 family) is likely the most realistic way to build a MuMax3-compatible “full quantum” module that materially outperforms purely classical LLB–Boltzmann on selected cluster-level problems. The extra qubits, improved coherence, and built-in error-mitigation workflows in the IBM Runtime stack give more headroom for larger local spin–bath registers and slightly deeper QITE circuits than on earlier ~133-qubit Heron-generation devices, while still staying within NISQ-era limits when combined with zero-noise extrapolation, measurement-error mitigation, and related QEM techniques.

* **Conceptually best:**
  SSE-, QMS-, and TFD-based methods can, in principle, still deliver higher physical fidelity and broader applicability (especially for strongly correlated spin models, non-Markovian baths, and finite-temperature entanglement and correlation functions) than any Ehrenfest-style QITE solver. However, even on a 156-qubit Heron-class processor with QEM, they remain more demanding in qubits, circuit depth, sampling overhead, and variational optimization complexity, so they are better viewed as high-accuracy cluster-level or equilibrium-only subroutines that complement, rather than replace, an Ehrenfest–LLB–Boltzmann+QITE+QEM “MuMax3-quantum” module in the near term.

---

To predict Curie temperature for magnetic materials:
- Solve the stochastic Landau Lifshitz Gilbert (LLG) or the Landau Lifshitz Bloch (LLB) equations. 
  - Mathematical description is available in the latter half of this markdown. 
- MuMax3 solves a micromagnetic effective field model whose equation of motion is the (stochastic) Landau–Lifshitz–Gilbert (LLG) equation.
- The LLG equation itself is a phenomenological, coarse-grained model derived from underlying atomic spins and quantum mechanics, so the whole MuMax3 framework is “effective” by construction.
- The stochastic LLG equation can be solved with MuMax3 if done carefully.
- This can be performed in Python on Google Colab for example, especially on the GPU for working with larger data. 

Other notes: 
- Solving the standard LLG equation by itself cannot predict the Curie temperature. 
- MuMax adapted with Python (optionally in Google Colab) can be used to solve the stochastic LLG for predicting Curie temperature and other temperature dependent dynamics, but it cannot solve the LLB equation.
- To solve the LLB equation for magnetic modeling can be performed in Python, and if written carefully, can also mapped to quantum gates in a library such as Qiskit to run on a backend of interest. 

---

## Typical high performance approaches for solving magnetism dynamics and other considerations 

```
Magnetization Dynamics (micromagnetic scale)
└─ MuMax3
   ├─ Core engine (C/CUDA + Go)
   │   ├─ Solves LLG/ sLLG
   │   │   ├─ Deterministic LLG (T = 0)
   │   │   └─ Stochastic LLG (Temp > 0, Brown thermal field)
   │   ├─ Effective fields
   │   │   ├─ Exchange, demag, anisotropy, DMI
   │   │   ├─ Spin-transfer torque (STT), spin–orbit torque (SOT)
   │   │   └─ Thermal field H_therm(T, ΔV, Δt)
   │   └─ Time-steppers
   │       ├─ Fixed-step schemes
   │       └─ Adaptive sLLG schemes (Leliaert et al.)
   ├─ Internal scripting (.mx3, Go-like)
   │   ├─ Set mesh, regions, Msat, Aex, Ku, …
   │   ├─ Set Temp, ThermSeed, B_ext(t), J(t)
   │   └─ Run/ Minimize/ Relax/ Save
   └─ External ecosystem
       ├─ Python wrappers (e.g., mumax3c, custom scripts)
       │   ├─ Generate .mx3 scripts programmatically
       │   ├─ Launch MuMax3 jobs (local or cluster)
       │   └─ Parse OVF/ table output, compute observables
       ├─ Other tools (ParaView, OVF viewers)
       └─ Multiscale coupling (atomistic spin, LLB, Qiskit, etc.)
```

## Ideas for a kind of "quantum MuMax3" simulation workflow 

```
Goal: “MuMax3 full quantum solver” for magnetization dynamics
    (backed by a 156-qubit IBM Heron r2 processor with quantum error mitigation)

├─ 0. Classical baseline
│  └─ LLB–Boltzmann
│     ├─ LLB for macrospin m(t) (transverse + longitudinal damping)
│     └─ Boltzmann for electrons/phonons (ultrafast kinetics, spin-flip channels)
│
├─ 1. Quantum-augmented effective models
│  ├─ qLLB (Nieves et al.) → quantum corrections to LLB rates/statistics
│  ├─ Quantum LLG/ q-LLG (Wieser; Liu; Azimi-Mousolou, etc.)
│  └─ Ehrenfest coupling to classical fields
│     (Ehrenfest–LLB–Boltzmann: quantum spin cluster + classical micromagnetic field)
│
└─ 2. NISQ-compatible quantum solvers
      (156-qubit IBM Heron r2 + QEM layer via Qiskit Runtime)

   ├─ Hardware + QEM envelope
   │   ├─ Processor: IBM Heron r2, 156-qubit heavy-hex lattice,
   │   │   tunable couplers, TLS-mitigation-enhanced coherence (e.g. ibm_fez). 
   │   └─ Quantum error mitigation layer:
   │       • Noise characterization on Heron (sparse Pauli-Lindblad models, etc.)
   │       • Techniques: Pauli twirling/ randomized compiling, ZNE, PEC,
   │         measurement error mitigation, hybrid QEM protocols.
   │       • Implemented through Qiskit Runtime error-mitigation passes.
   │
   ├─ 2A. QITE-based Ehrenfest–LLB–Boltzmann (primary MuMax3-quantum module)
   │   ├─ Encode local Ehrenfest–qLLB–Boltzmann generator as a linear operator
   │   │   on a spin + bath register (subset of the 156 qubits).
   │   ├─ Use QITE-style circuits to prepare steady/ thermal states and
   │   │   short real-time steps for ⟨m(t)⟩, ⟨E(t)⟩.
   │   └─ Run on Heron r2 with QEM:
   │       • Compile to heavy-hex topology, optimize depth.
   │       • Apply ZNE + Pauli twirling and/or PEC to magnetization observables.
   │
   ├─ 2B. SSE-based trajectories on Heron r2
   │   ├─ Variational SSE (Endo; Lan): simulate stochastic Schrödinger trajectories
   │   │   for local quantum spin clusters (tens of spins + ancillas).
   │   ├─ Optionally embed Generalized Quantum Master Equation (GQME) kernels
   │   │   to capture non-Markovian baths.
   │   └─ Use QEM on expectation values from trajectory ensembles
   │       (ZNE/PEC + measurement mitigation).
   │
   ├─ 2C. QMS (quantum Metropolis) on Heron r2
   │   ├─ Quantum-Metropolis(-Hastings) sampling for Gibbs states of spin
   │   │   Hamiltonians relevant to micromagnetics (Heisenberg + anisotropy, etc.).
   │   ├─ Use the 156-qubit space to host moderate-size quantum spin lattices
   │   │   whose thermal states are hard for classical Monte Carlo (sign problems).
   │   └─ Apply QEM to thermal expectation values (magnetization curves, χ(T), etc.).
   │
   └─ 2D. TFD-based finite-T simulation on Heron r2
       ├─ Variational TFD preparation (Sagastizabal, Su, Zhu style circuits)
       │   for embedded spin clusters, using up to ~1/2 of the 156 qubits per copy.
       ├─ Real-time evolution of TFD states to compute finite-temperature
       │   correlation functions and dynamical susceptibilities.
       └─ QEM applied to TFD preparation and time evolution (twirling, ZNE, PEC),
           focusing on thermodynamic and linear-response observables.

└─ 3. Future fault-tolerant/ post-NISQ directions
   ├─ Larger qLLG/ qLLB clusters with exact SSE or full master-equation solvers
   │  (multiple Heron-class nodes or error-corrected logical qubits).
   ├─ Full QMS phase diagrams for realistic 2D/3D spin Hamiltonians
   │  (true quantum criticality and Curie-temperature benchmarks).
   └─ Hierarchical embedding:
      ├─ Quantum clusters (solved on logical qubits) providing qLLB/ kinetic kernels.
      └─ These kernels drive a classical MuMax3-like grid at device scale.
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
| **Kinetic‑aware Ehrenfest–LLB–Boltzmann (E‑LLB‑B)** | **Yes** ($\Gamma_\parallel$($T$)) | **Yes** with detailed-balance Lindbladian anchoring | **Yes** (full $P(m)$ dynamics enables $\chi$ and Binder) | **Yes** via collision-integral-based rates (qLLB-style) | **Strong** (robust $T_c$ from $m(T)$, $\chi(T)$, Binder $U_4$) | Broad: equilibrium + kinetics + finite-size scaling |
| **Ehrenfest–LLB (GKSL, expectation-level)** | **Yes** (via Lindblad rates) | **Yes** (KMS/Davies choice) | Moments only; sample trajectories for chi and Binder | Indirect (time-dependent rates) | **Good** if rates are calibrated from microscopic models | Medium–high: efficient but does not evolve full $P(m)$ |

| Approach                                     | Physical content                                                                                                                 | Strengths for a “MuMax3-quantum” solver                                                                                                                                                                                                                                                                                   | Main limitations on 156-qubit Heron NISQ (r2)                                                                                                                                                                                                                                                                                             | How it compares to pure LLB-Boltzmann                                                                                                                                                                                         |
| -------------------------------------------- | -------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| **Classical LLB-Boltzmann**                  | Classical macrospin + Boltzmann kinetics; quantum physics only in parameters                                                     | Mature, efficient; already used for ultrafast magnetization and near-($T_C$) physics; easy to scale to full MuMax3-scale grids                                                                                                                                                                                              | No direct use of Heron; limited purely by classical resources. Misses quantum entanglement and quantum statistics of spins; macrospin approximation only                                                                                                                                                                                  | Baseline                                                                                                                                                                                                                      |
| **Ehrenfest-LLB-Boltzmann + QITE + QEM**     | qLLB-derived macrospin coupled to Boltzmann baths; solved via QITE on a quantum register; Ehrenfest coupling to classical fields | Brings in quantum statistics and local spin–bath correlations; leverages QITE as a PDE/ steady-state solver; QEM + dissipation can stabilize deeper circuits on Heron r2; 156 qubits and improved coherence/TLS mitigation allow larger local clusters and somewhat deeper circuits than 133-qubit Heron r1  | Still uses an Ehrenfest mean-field treatment; limited by NISQ depth, ~($10^{-3}$)–($10^{-2}$) 2-qubit error rates and heavy-hex connectivity (SWAP overhead); QEM (ZNE/PEC/readout mitigation) adds sampling and classical post-processing overhead; quantum speedup not guaranteed yet on 156-qubit devices                | **Clear upgrade**: more quantum-faithful magnetization dynamics and potentially better scaling for local quantum clusters, now with more qubits and better coherence headroom than on 133-qubit Heron r1                      |
| **SSE-based quantum trajectories (VQS-SSE)** | Full quantum trajectories for spins with dissipators; master-equation level                                                      | Captures fluctuations, jumps, non-Markovian effects; directly implements qLLG/ qLLB-like open dynamics; 156 qubits allow somewhat larger clusters (more spins + ancillas) per trajectory than on 133-qubit hardware                                                                                                      | Circuits are deeper and stochastic; even with Heron r2’s improved fidelities and TLS mitigation, depth + trajectory sampling burden remain large; embedding many such clusters into a full micromagnetic grid is still challenging on a single 156-qubit NISQ device                                                                      | **More accurate** at the local quantum-cluster level but likely **still less practical** on current 156-qubit Heron systems for large-scale MuMax3-style grids                                                                |
| **Quantum-Metropolis/ QMS**                 | True Gibbs states of quantum spin Hamiltonians; avoids sign problem                                                              | Excellent for equilibrium magnetism, phase diagrams, accurate ($T_C$) in strongly correlated regimes; 156 qubits can host larger spin lattices per sample than 133-qubit chips, improving finite-size scaling                                                                                                               | Requires phase-estimation-like subroutines and deep circuits that remain difficult even on Heron r2; QEM overhead further increases sampling cost; better suited to off-line equilibrium studies on modest-size clusters rather than continuous micromagnetic time-stepping                                                               | **Better equilibrium physics** than LLB-Boltzmann, especially for strongly correlated spin models, but not directly a time-evolution engine and still resource-heavy on 156-qubit NISQ hardware                               |
| **TFD-based VQA**                            | Purified thermal states + unitary real-time evolution                                                                            | Rigorous thermal correlators; natural way to study finite-($T$) quantum spin dynamics and entanglement; 156 qubits allow TFD of embedded spin clusters (e.g. ~60–70 physical spins + ancillas) with improved coherence and QEM support                                                                                      | Doubles system qubits for the TFD copy, so usable cluster size is effectively halved; variational TFD preparation and time-evolution circuits are deep and sensitive to noise; although Heron r2’s TLS mitigation and QEM help, current demonstrations are on smaller systems and scaling to the full 156-qubit space remains challenging | **Most accurate finite-(T) description** in principle (for the cluster it fits), but expensive; ideal as a high-fidelity cluster-level module feeding effective parameters back into a larger LLB-Boltzmann/ MuMax3 workflow |

---

## GPU calculated result for Curie temperature prediction 

<img width="1864" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/23c695e1-4b90-4a50-bbfd-f44ceb211da7" />


## Real QPU calculated result (without error mitigation) for Curie temperature prediction 

<img width="1864" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/b3f3bb24-bcda-4a7a-a094-70f9ddf6f22a" />

---

## Relevant equations:

### Notation

* Magnetization (unit): ( $\boldsymbol{m}=\boldsymbol{M}/M_s$ ), ( $m=\lvert\boldsymbol{m}\rvert$ ).
* Effective field: ( $\boldsymbol{B}_{\mathrm{eff}}=-\partial \mathcal{F}/\partial \boldsymbol{m}$ ), with Zeeman, exchange, anisotropy, Dzyaloshinskii–Moriya (DMI), etc.
* Quantum expectation (site (i)): ( $\boldsymbol{s}_i(t)=\mathrm{Tr}!\big[\rho(t),\boldsymbol{\sigma}_i\big]$ ).
* Unit field direction: ( $\hat{\boldsymbol{b}}=\boldsymbol{B}*{\mathrm{eff}}/\lvert \boldsymbol{B}*{\mathrm{eff}}\rvert$ ).
* Rates: longitudinal ( $\Gamma_\parallel=1/T_1$ ), transverse ( $\Gamma_\perp=1/T_2$ ); LLB coefficients ( $\alpha_\parallel(T),\alpha_\perp(T)$ ).
* Equilibrium magnetization: ( $m_{\mathrm{eq}}(T)$ ) (or ( $m_{\mathrm{eq}}(T,D)$ ) if a finite-size parameter ( $D$ ) is included).
* Pauli vector: ( $\boldsymbol{\sigma}_i=(\sigma_i^x,\sigma_i^y,\sigma_i^z)$ ).
* Density matrix: ( $\rho(t)$ ); pure state trajectories: ( $\lvert\psi(t)\rangle$ ).
* Bohr frequencies: ( $\omega_{mn}=(E_m-E_n)/\hbar$ ), inverse temperature ( $\beta=1/k_{\mathrm B}T$ ).

---

## 0) Ehrenfest theorems (quantum → classical bridge)

### 0.1 Canonical Ehrenfest theorem (position–momentum pair)

For a particle of mass ( $m$ ) with Hamiltonian ( $H=\hat{p}^2/(2m)+V(\hat{x},t)$ ):

$$
\frac{d}{dt}\langle \hat{x}\rangle
= \frac{1}{m},\langle \hat{p}\rangle,
$$

$$
\frac{d}{dt}\langle \hat{p}\rangle
= -\left\langle \frac{\partial V(\hat{x},t)}{\partial x}\right\rangle
= \langle \hat{F}\rangle.
$$

These are the “canonical” Ehrenfest relations: the expectation values obey Newton-like equations whenever the wave packet remains sufficiently localized. ([Physical Review Links][1])

### 0.2 General Ehrenfest theorem (any operator ( $\hat{A}$ ))

For a general (possibly time-dependent) operator ( $\hat{A}(t)$ ):

$$
\frac{d}{dt}\langle \hat{A}\rangle
= \left\langle\frac{\partial \hat{A}}{\partial t}\right\rangle

* \frac{1}{i\hbar},\big\langle[\hat{A},\hat{H}]\big\rangle.
  $$

This is the **general Ehrenfest theorem**, valid in Schrödinger or Heisenberg pictures. Choosing appropriate ( $\hat{A}$ ) (e.g., spin components, magnetization operators) yields all the Ehrenfest magnetization equations below. ([Physical Review Links][1])

---

## 1) Landau–Lifshitz / Gilbert and stochastic LLG

### 1.1 Classical original (Gilbert and LL forms)

**Gilbert form**:

$$
\frac{d\boldsymbol{m}}{dt}
= -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}

* \frac{\alpha}{m},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt}.
  $$

**Landau–Lifshitz (LL) form**:

$$
\frac{d\boldsymbol{m}}{dt}
= -\gamma',\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}

* \gamma'\alpha,\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}\big),
  \qquad
  \gamma'=\frac{\gamma}{1+\alpha^2}.
  $$

### 1.2 Stochastic LLG (thermal field)

$$
\frac{d\boldsymbol{m}}{dt}
= -\gamma,\boldsymbol{m}\times\big(\boldsymbol{B}*{\mathrm{eff}}+\boldsymbol{B}*{\mathrm{th}}(t)\big)

* \frac{\alpha}{m},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt},
  $$

$$
\big\langle B_{\mathrm{th},\mu}(t),B_{\mathrm{th},\nu}(t')\big\rangle
= 2D,\delta_{\mu\nu},\delta(t-t'),
\qquad
D\propto\frac{\alpha k_{\mathrm B}T}{\gamma M_s V}.
$$

This is the **stochastic LLG** underlying MuMax3-style micromagnetics.

### 1.3 Ehrenfest quantum Landau–Lifshitz (E-qLL)

Take a quantum spin (or spin cluster) with Hamiltonian ( $\hat{H}$ ) and define
( $\boldsymbol{s}(t)=\tfrac12\langle\boldsymbol{\sigma}(t)\rangle$ ). From the general Ehrenfest theorem:

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{1}{i\hbar},\big\langle[\tfrac12\boldsymbol{\sigma},\hat{H}]\big\rangle
\equiv \boldsymbol{s}\times\boldsymbol{\Omega}_{\mathrm{eff}},
$$

where ( $\boldsymbol{\Omega}*{\mathrm{eff}}$ ) is an operator-averaged precession vector (for ( $\hat{H}= -\tfrac{\hbar\gamma}{2}\boldsymbol{B}*{\mathrm{eff}}!\cdot!\boldsymbol{\sigma}$ one recovers the LL precession term).

Adding phenomenological transverse damping at the expectation level gives

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}}
-\Gamma*{\perp}\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big),
\qquad
\hat{\boldsymbol{b}} = \frac{\boldsymbol{B}*{\mathrm{eff}}}{\lvert\boldsymbol{B}*{\mathrm{eff}}\rvert},
$$

which is the **Ehrenfest quantum Landau–Lifshitz (E-qLL)** equation (transverse damping only).

### 1.4 Ehrenfest quantum Landau–Lifshitz–Gilbert (E-qLLG)

A **quantum analog of LLG (q-LLG)** can be formulated for the density operator ( $\rho$ ) of a spin system (Liu *et al.*). ([Physical Review Links][2])

At the level of ( $\boldsymbol{s}(t)=\mathrm{Tr}[\rho(t),\boldsymbol{\sigma}/2]$ ), one obtains an Ehrenfest-type qLLG:

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}

* \alpha_{\mathrm{q}};\boldsymbol{s}\times\frac{d\boldsymbol{s}}{dt},
  $$

or, equivalently in LL form,

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}

* \frac{2\alpha_{\mathrm{q}}}{\hbar},
  \boldsymbol{s}\times\big(\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}\big),
  $$

with a **quantum-derived damping** coefficient ( $\alpha_{\mathrm{q}}$ ) constrained by the underlying q-LLG master equation. ([Physical Review Links][2])

This is the **Ehrenfest quantum Landau–Lifshitz–Gilbert (E-qLLG)** equation.

### 1.5 Ehrenfest quantum Landau–Lifshitz–Gilbert–Slonczewski (E-qLLGS)

Include a **Slonczewski spin-transfer torque (STT)** term with polarization ( $\hat{\boldsymbol{p}}$ ):

$$
\left.\frac{d\boldsymbol{s}}{dt}\right|_{\mathrm{STT}}
= -\gamma a_J,\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big),
$$

so that the full Ehrenfest **qLLG–Slonczewski (E-qLLGS)** equation becomes

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}

* \frac{2\alpha_{\mathrm{q}}}{\hbar},
  \boldsymbol{s}\times\big(\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}\big)
* \gamma a_J,\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big).
  $$

At the GKSL level this can be realized by adding spin-current-induced Lindblad channels or effective non-Hermitian terms coupled to a reservoir of polarized conduction electrons.

### 1.6 Adapted (Qiskit-ready, LL/LLG sector)

**Hamiltonian (single macrospin or per site (i))**

$$
\hat{H}(t)= -\frac{\hbar\gamma}{2},\boldsymbol{B}_{\mathrm{eff}}(t)\cdot\boldsymbol{\sigma}
\quad\text{(encode as a sum of Pauli terms).}
$$

**Optional LL/LLG-type damping on simulators**

Add appropriate Lindblad operators ( $L^\pm,L^z$ ) (see §4) so that the Ehrenfest equations reproduce E-qLL or E-qLLG. This is directly compatible with a Lindblad-based simulator or a QITE / VQS approximation to Lindblad dynamics in Qiskit ≥ 2.2. ([ADS][3])

---

## 2) Landau–Lifshitz–Bloch (LLB) and dynamic extensions

### 2.1 Classical original LLB (single-sublattice macrospin)

$$
\frac{d\boldsymbol{m}}{dt}
= -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}

* \frac{\gamma,\alpha_{\parallel}(T)}{m^2},\big(\boldsymbol{m}!\cdot!\boldsymbol{B}_{\mathrm{eff}}\big),\boldsymbol{m}

- \frac{\gamma,\alpha_{\perp}(T)}{m^2},
  \boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}\big).
  $$

A standard longitudinal field component inside ( $\boldsymbol{B}_{\mathrm{eff}}$ ) is

$$
\boldsymbol{B}*{\parallel}
= \frac{1}{\chi*{\parallel}(T)}
\left(1-\frac{m^2}{m_{\mathrm{eq}}^2(T)}\right)\boldsymbol{m},
\quad\Rightarrow\quad
\lvert \boldsymbol{m}\rvert \to m_{\mathrm{eq}}(T).
$$

### 2.2 Dynamic LLB (time-dependent fields and parameters)

Allow explicit time dependence (e.g. ultrafast laser pumping, evolving temperature ( $T(t)$ )):

$$
\frac{d\boldsymbol{m}}{dt}
= -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}(\boldsymbol{m},t)

* \frac{\gamma,\alpha_{\parallel}(T(t))}{m^2},\big(\boldsymbol{m}!\cdot!\boldsymbol{B}_{\mathrm{eff}}(\boldsymbol{m},t)\big),\boldsymbol{m}

- \frac{\gamma,\alpha_{\perp}(T(t))}{m^2},
  \boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}(\boldsymbol{m},t)\big).
  $$

This is the **dynamic Landau–Lifshitz–Bloch (d-LLB)** equation.

### 2.3 Dynamic LLB–Slonczewski (d-LLBS)

Add Slonczewski spin-transfer torque to d-LLB:

$$
\left.\frac{d\boldsymbol{m}}{dt}\right|_{\mathrm{STT}}
= -\gamma a_J(t),\boldsymbol{m}\times\big(\boldsymbol{m}\times\hat{\boldsymbol{p}}\big),
$$

so

$$
\frac{d\boldsymbol{m}}{dt}
= \text{d-LLB drift from §2.2}

* \gamma a_J(t),\boldsymbol{m}\times\big(\boldsymbol{m}\times\hat{\boldsymbol{p}}\big).
  $$

This is the **dynamic Landau–Lifshitz–Bloch–Slonczewski (d-LLBS)** equation.

### 2.4 Ehrenfest classical LLB (deterministic mean)

Replace ( $\boldsymbol{m}$ ) by quantum spin expectation ( $\boldsymbol{s}$ ) and write the LLB in Ehrenfest form:

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}

* \Gamma_{\perp}(T)\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big)
* \Gamma_{\parallel}(T)\Big[(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}},
  $$

with ( $\hat{\boldsymbol{b}}=\boldsymbol{B}*{\mathrm{eff}}/\lvert\boldsymbol{B}*{\mathrm{eff}}\rvert$ ) and classical LLB rates related to ( $\alpha_{\parallel,\perp}(T)$ ).

### 2.5 Quantum Landau–Lifshitz–Bloch (qLLB) — density-matrix form

Nieves *et al.* derive a **quantum LLB (qLLB)** equation from a spin–phonon / spin–electron density-matrix treatment. ([arXiv][4])

At the master-equation level:

$$
\dot{\rho}
= -\frac{i}{\hbar}[\hat{H},\rho]

* \sum_{\mu}\mathcal{D}[L_\mu]\rho,
  $$

with spin-flip and dephasing jump operators (e.g., ( $L^\pm\propto\sigma^\pm$ ), ( $L^z\propto\sigma^z$ )) chosen so that the Ehrenfest equation for ( $\boldsymbol{s}(t)=\mathrm{Tr}[\rho,\boldsymbol{\sigma}/2]$ ) reproduces an LLB-type structure, but with **quantum-corrected rates**:

$$
\alpha_{\parallel}(T)
= \lambda,\frac{2T}{3T_C},\frac{2q_s}{\sinh(2q_s)},
\qquad
\alpha_{\perp}(T)
= \lambda\left[\frac{\tanh q_s}{q_s}-\frac{2T}{3T_C}\right],
\qquad
q_s=\frac{\mu H_{\mathrm{MFA}}}{k_{\mathrm B}T}.
$$

### 2.6 Ehrenfest quantum LLB (E-qLLB)

From qLLB, the Ehrenfest equation for ( $\boldsymbol{s}$ ) yields **Ehrenfest quantum LLB**:

$$
\frac{d\boldsymbol{s}}{dt}
= \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}

* \Gamma_{\perp}^{\mathrm{q}}(T)\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big)
* \Gamma_{\parallel}^{\mathrm{q}}(T)\Big[(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}^{\mathrm{q}}(T)\Big]\hat{\boldsymbol{b}},
  $$

with ( $\Gamma_{\parallel,\perp}^{\mathrm{q}}(T)$ ) and ( $m_{\mathrm{eq}}^{\mathrm{q}}(T)$ ) obtained from qLLB.

### 2.7 Quantum dynamic LLB (q-dLLB) and E-q-dLLB

Allow explicit time dependence ( $T(t)$, ( $\boldsymbol{B}_{\mathrm{eff}}(t)$ ), etc.) in qLLB:

$$
\dot{\rho}(t)
= -\frac{i}{\hbar}[\hat{H}(t),\rho(t)]

* \sum_{\mu}\mathcal{D}[L_\mu(t)]\rho(t)
  \quad\Rightarrow\quad
  \frac{d\boldsymbol{s}}{dt}
  = \text{E-qLLB drift with }T\to T(t),;\boldsymbol{B}*{\mathrm{eff}}\to\boldsymbol{B}*{\mathrm{eff}}(t),
  $$

giving **quantum dynamic LLB (q-dLLB)** for ( $\rho$ ) and **Ehrenfest quantum dynamic LLB (E-q-dLLB)** for ( $\boldsymbol{s}$ ).

### 2.8 Quantum dynamic LLB–Slonczewski (q-dLLBS) and E-q-dLLBS

Add STT-type channels / Hamiltonian terms to q-dLLB:

$$
\dot{\rho}(t)
= -\frac{i}{\hbar}\big[\hat{H}(t)+\hat{H}_{\mathrm{STT}}(t),\rho(t)\big]

* \sum_{\mu}\mathcal{D}[L_\mu(t)]\rho(t),
  $$

with ( $\hat{H}_{\mathrm{STT}}\propto a_J(t),\hat{\boldsymbol{S}}\cdot(\hat{\boldsymbol{S}}\times\hat{\boldsymbol{p}})$ ) in an appropriate spin representation. At the Ehrenfest level this yields **Ehrenfest quantum dynamic LLB–Slonczewski (E-q-dLLBS)**:

$$
\frac{d\boldsymbol{s}}{dt}
= \text{E-q-dLLB drift of §2.7}

* \gamma a_J(t),\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big).
  $$

---

## 3) Kinetic-aware Ehrenfest–Boltzmann models on the Bloch ball

Introduce a distribution ( $f(\boldsymbol{m},t)$ ) (or ( $f(\boldsymbol{s},t)$ )) on the Bloch ball and write a drift–diffusion–collision kinetic equation:

$$
\partial_t f

* \nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{A}(\boldsymbol{m},T,t),f\big]
  = \nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{D}(\boldsymbol{m},T,t),\nabla_{\boldsymbol{m}} f\big]
* \mathcal{C}[f],
  $$

where ( $\mathcal{C}[f]$ ) encodes magnon/electron/phonon collisions and ( $\boldsymbol{D}$ ) satisfies fluctuation–dissipation.

### 3.1 Ehrenfest–LLB–Boltzmann (E-LLB-B)

Use **LLB drift**:

$$
\boldsymbol{A}*{\mathrm{LLB}}(\boldsymbol{m},T,t)
= \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}}

* \Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big)
* \Gamma_{\parallel}(T)\Big[(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}}.
  $$

With ( $\Gamma_{\parallel,\perp}$ ) and ( $m_{\mathrm{eq}}(T)$ ) obtained from **classical** LLB, this defines the **Ehrenfest–LLB–Boltzmann (E-LLB-B)** equation.

### 3.2 Ehrenfest–quantum-LLB–Boltzmann (E-qLLB-B)

Replace classical LLB rates by **qLLB rates** (§2.5):

$$
\boldsymbol{A}*{\mathrm{qLLB}}(\boldsymbol{m},T,t)
= \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}}

* \Gamma_{\perp}^{\mathrm{q}}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big)
* \Gamma_{\parallel}^{\mathrm{q}}(T)\Big[(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}^{\mathrm{q}}(T)\Big]\hat{\boldsymbol{b}}.
  $$

Then **Ehrenfest quantum LLB–Boltzmann (E-qLLB-B)** is the kinetic equation with ( $\boldsymbol{A}_{\mathrm{qLLB}}$ ).

### 3.3 Ehrenfest–LL–Boltzmann (E-LL-B) and E-qLL-B

Use **LL drift only**:

$$
\boldsymbol{A}*{\mathrm{LL}}(\boldsymbol{m},T,t)
= \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}}

* \Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big),
  $$

which lacks a longitudinal channel, so it cannot alone produce the amplitude collapse ( $m\to 0$ ) at ( $T_C$ ). With quantum-derived ( $\Gamma_{\perp}^{\mathrm{q}}(T)$ ) this becomes **Ehrenfest quantum LL–Boltzmann (E-qLL-B)**.

### 3.4 Dynamic and Slonczewski variants (d-E-LLB-B, d-E-qLLB-B-S)

For the **dynamic** and **Slonczewski** cases:

* **Dynamic Ehrenfest–LLB–Boltzmann (d-E-LLB-B)**: use ( $\boldsymbol{A}_{\mathrm{LLB}}(\boldsymbol{m},T(t),t)$ ).
* **Dynamic Ehrenfest–quantum-LLB–Boltzmann (d-E-qLLB-B)**: use ( $\boldsymbol{A}_{\mathrm{qLLB}}(\boldsymbol{m},T(t),t)$ ).
* **Dynamic Ehrenfest–quantum-LLB–Boltzmann–Slonczewski (d-E-qLLB-B-S)**: add STT drift
  ( $-\gamma a_J(t),\boldsymbol{m}\times(\boldsymbol{m}\times\hat{\boldsymbol{p}})$ ) to ( $\boldsymbol{A}_{\mathrm{qLLB}}$ ).

These are your **Ehrenfest quantum Landau–Lifshitz–Bloch–Boltzmann**, **Ehrenfest quantum Landau–Lifshitz–Boltzmann**, **dynamic** and **Slonczewski-extended** Boltzmann couplings.

---

## 4) Lindblad (GKSL) backbone and thermal anchoring (quantum-native)

### 4.1 GKSL master equation

For all quantum-native versions (qLL, qLLG, qLLB, q-dLLB, q-dLLBS), the density matrix obeys a GKSL equation:

$$
\dot{\rho}
= -\frac{i}{\hbar}[\hat{H},\rho]

* \sum_\mu\Big(L_\mu,\rho,L_\mu^\dagger
  -\tfrac12{L_\mu^\dagger L_\mu,\rho}\Big)
  \equiv -\frac{i}{\hbar}[\hat{H},\rho]
* \sum_\mu \mathcal{D}[L_\mu]\rho.
  $$

With suitable ( $L_\mu$ ), **Landau–Lifshitz damping can be derived directly from Lindbladian dissipation**, establishing LL/LLG/LLB as Ehrenfest limits of GKSL dynamics. ([Physical Review Links][2])

### 4.2 Detailed balance (KMS) and ( $T_1/T_2$ ) dictionary

For a two-level splitting ( $\hbar\omega$ ) and spin-flip Lindblad operators:

$$
L^{-} = \sqrt{\gamma_{\downarrow}},\sigma^{-},
\quad
L^{+} = \sqrt{\gamma_{\uparrow}},\sigma^{+},
\quad
L^{z} = \sqrt{\gamma_{\phi}},\sigma^{z},
$$

imposing **KMS detailed balance**:

$$
\frac{\gamma_{\uparrow}}{\gamma_{\downarrow}}
= e^{-\beta\hbar\omega},
\qquad
\Gamma_{\parallel} = \gamma_{\downarrow}+\gamma_{\uparrow},
\qquad
\Gamma_{\perp} = \frac{\Gamma_{\parallel}}{2}+\gamma_{\phi},
$$

gives the thermal ( $T_1/T_2$ ) dictionary underlying LLB/qLLB.

### 4.3 Quantum-native mappings for QITE, SSE, and TFD

All the GKSL-level equations in §§1–3 can be mapped to **quantum-circuit-compatible algorithms**:

1. **QITE / VQS for Lindblad**
   Treat GKSL as a first-order linear differential equation and apply **Variational Quantum Simulation (VQS)** or **Quantum Imaginary-Time Evolution (QITE)** extensions to open systems. ([ADS][3])
   *Implementation:*

   * Pauli-decompose ( $\hat{H}$ ) and ( $L_\mu$ ) (cf. §5).
   * Use Qiskit ≥ 2.2’s primitive-based VQS/QITE infrastructure with Quantum Error Mitigation (ZNE, MEM, etc.) on a 156-qubit Heron backend.

2. **Stochastic Schrödinger Equation (SSE) trajectories**
   Unravel GKSL into a **stochastic Schrödinger equation**:

   $$
   d\lvert\psi(t)\rangle
   = -\frac{i}{\hbar}\hat{H}_{\mathrm{eff}}(t)\lvert\psi(t)\rangle,dt

   * \sum_\mu \Big(\frac{L_\mu}{\sqrt{p_\mu(t)}}-\mathbb{I}\Big)\lvert\psi(t)\rangle,dN_\mu(t),
     $$

   where ( $\hat{H}*{\mathrm{eff}}$ ) includes non-Hermitian contributions, and ( $dN*\mu(t)$ ) are Poisson increments. Variational SSE implementations on NISQ devices have been developed for time-local master equations and non-Markovian baths. ([ADS][3])

3. **Thermofield-Double (TFD)-based VQAs**
   For finite-temperature qLLB/q-dLLB/q-dLLBS, one can represent thermal states as **TFD states** on a doubled Hilbert space and optimize a variational circuit to approximate ( $\lvert\mathrm{TFD}(\beta)\rangle$ ), then evolve unitarily. ([Physical Review Links][1])

   *Implementation:*

   * Prepare TFD using QAOA-style or VQE-style ansatz.
   * Evolve with qLL/qLLG/qLLB Hamiltonians (plus ancilla-encoded dissipation).
   * Use QEM (ZNE, MEM, possibly probabilistic error cancellation on small subsystems) on a 156-qubit Heron processor.

In all cases, the Ehrenfest magnetization and kinetic equations in §§1–3 are recovered by taking expectation values of the quantum trajectories or density matrix, thus giving you quantum-native realizations of:

* canonical and general Ehrenfest equations,
* E-qLL, E-qLLG, E-qLLGS,
* E-qLLB, E-q-dLLB, E-q-dLLBS,
* and their Boltzmann-coupled variants.

---

## 5) Adapted Hamiltonian building blocks (for Qiskit encodings)

Nearest-neighbor Heisenberg + Zeeman + DMI on a qubit graph:

$$
\hat{H}
= -\frac{\hbar\gamma}{2}\sum_i \boldsymbol{B}_{\mathrm{eff},i}!\cdot!\boldsymbol{\sigma}_i

* \sum_{\langle i,j\rangle}\big(J_x X_iX_j+J_y Y_iY_j+J_z Z_iZ_j\big)

- \sum_{\langle i,j\rangle}\boldsymbol{D}_{ij}!\cdot!\big(\boldsymbol{\sigma}_i\times\boldsymbol{\sigma}_j\big)
- \hat{H}*{\mathrm{ani}} + \hat{H}*{\mathrm{STT}}.
  $$

**Lindblad channels (per site (i))**

$$
L_i^{-}=\sqrt{\gamma_{\downarrow,i}},\sigma_i^{-},
\quad
L_i^{+}=\sqrt{\gamma_{\uparrow,i}},\sigma_i^{+},
\quad
L_i^{z}=\sqrt{\gamma_{\phi,i}},\sigma_i^{z},
$$

with KMS-consistent rates as in §4.2. These primitives reproduce **Ehrenfest LL/LLG/LLB/qLLB** drifts when inserted into the GKSL backbone (§4.1), and their Pauli decompositions are directly usable in Qiskit circuits.

---

## 6) Index of equations (chronological / conceptual)

* **Canonical Ehrenfest theorem:** §0.1.

* **General Ehrenfest theorem:** §0.2.

* **LLG / LL (classical):** §1.1.

* **sLLG:** thermal-field form + white-noise correlator, §1.2.

* **Ehrenfest quantum LL (E-qLL):** §1.3.

* **Ehrenfest quantum LLG (E-qLLG):** §1.4.

* **Ehrenfest quantum LLG–Slonczewski (E-qLLGS):** §1.5.

* **Classical LLB:** §2.1.

* **Dynamic LLB (d-LLB):** §2.2.

* **Dynamic LLB–Slonczewski (d-LLBS):** §2.3.

* **Ehrenfest classical LLB:** §2.4.

* **Quantum LLB (qLLB, GKSL form):** §2.5.

* **Ehrenfest quantum LLB (E-qLLB):** §2.6.

* **Quantum dynamic LLB (q-dLLB) + E-q-dLLB:** §2.7.

* **Quantum dynamic LLB–Slonczewski (q-dLLBS) + E-q-dLLBS:** §2.8.

* **Ehrenfest–LLB–Boltzmann (E-LLB-B):** §3.1.

* **Ehrenfest–quantum-LLB–Boltzmann (E-qLLB-B):** §3.2.

* **Ehrenfest–LL–Boltzmann (E-LL-B) and E-qLL-B:** §3.3.

* **Dynamic Ehrenfest–(q)LLB–Boltzmann(+Slonczewski):** §3.4.

* **GKSL master equation (quantum-native backbone):** §4.1.

* **KMS + ( $T_1/T_2$ ) dictionary:** §4.2.

* **Quantum-native QITE / SSE / TFD mappings (Qiskit-compatible):** §4.3.

* **Adapted Hamiltonian + dissipators (Qiskit Pauli form):** §5.

---

### Symbol mini‑glossary

( $\gamma$ ): gyromagnetic ratio; ( $\alpha$ ): Gilbert damping; ( $\alpha_{\parallel,\perp}(T)$ ): LLB coefficients; ( $\Gamma_{\parallel,\perp}$ ): Ehrenfest rates; ( $m_{\mathrm{eq}}(T)$ ): equilibrium magnetization; ( $\chi_{\parallel}(T)$ ): longitudinal susceptibility; ( $\boldsymbol{\sigma}$ = $(X,Y,Z)$ ): Pauli vector; ( $J_{x,y,z}$ ): exchange couplings; ( $\boldsymbol{D} * {ij}$ ): DMI vector; ( $\gamma*{\uparrow,\downarrow,\phi}$ ): Lindblad rate parameters; ( $\beta=(k_{\mathrm B}T)^{-1}$ ).

---

References

1.  Atxitia, U., Hinzke, D. and Nowak, U. (2017) 'Fundamentals and applications of the Landau–Lifshitz–Bloch equation', *Journal of Physics D: Applied Physics*, 50(3), p. 033003. Available at: <https://doi.org/10.1088/1361-6463/50/3/033003>.
2.  Castin, Y., Dalibard, J. and Mølmer, K. (2008) *A Wave Function approach to dissipative processes*. arXiv:0805.4002. Available at: <https://doi.org/10.48550/arXiv.0805.4002>.
3.  Chen, C-F., Kastoryano, M., Brandão, F.G.S.L. and Gilyén, A. (2025) 'Efficient quantum thermal simulation', *Nature*, 646(8085), pp. 561–566. Available at: <https://doi.org/10.1038/s41586-025-09583-x>.
4.  Dalibard, J., Castin, Y. and Mølmer, K. (1992) 'Wave-function approach to dissipative processes in quantum optics', *Physical Review Letters*, 68(5), pp. 580–583. Available at: <https://doi.org/10.1103/PhysRevLett.68.580>.
5.  Devyaterikov, D.I., Proglyado, V.V., Zhaketov, V.D., Nikitenko, Y.V., Kondrat'ev, O.A., Pashaev, E.M., Subbotin, I.A., Zverev, V.I., Kravtsov, E.A. and Ustinov, V.V. (2021) 'Influence of Dimensional Effects on the Curie Temperature of Dy and Ho Thin Films', *Physics of Metals and Metallography*, 122(5), pp. 465–471. Available at: <https://doi.org/10.1134/S0031918X21050033>.
6.  Evans, R.F.L., Hinzke, D., Atxitia, U., Nowak, U., Chantrell, R.W. and Chubykalo-Fesenko, O. (2012) 'Stochastic form of the Landau-Lifshitz-Bloch equation', *Physical Review B*, 85(1), p. 014433. Available at: <https://doi.org/10.1103/PhysRevB.85.014433>.
7.  Garanin, D.A. (1998) *Fokker-Planck and Landau-Lifshitz-Bloch equations for classical ferromagnets*. arXiv:cond-mat/9805054. Available at: <https://doi.org/10.48550/arXiv.cond-mat/9805054>.
8.  Gokhale, S. and Manna, U. (2023) *Optimal control of the stochastic Landau-Lifshitz-Bloch equation*. arXiv:2305.10861. Available at: <https://doi.org/10.48550/arXiv.2305.10861>.
9.  Liang, J-M., Lv, Q-Q., Wang, Z-X. and Fei, S-M. (2023) 'Assisted quantum simulation of open quantum systems', iScience, 26(4), p. 106306. Available at: <https://doi.org/10.1016/j.isci.2023.106306>.
10. Liu, S.H., Behrendt, D.R., Legvold, S. and Good, R.H., Jr. (1959) 'Interpretation of Magnetic Properties of Dysprosium', *Physical Review*, 116(6), pp. 1464–1468. Available at: <https://doi.org/10.1103/PhysRev.116.1464>.
11. Meil, D., Evelt, M., Pfau, B., Kläui, M., Atxitia, U. and Nowak, U. (2020) *Thermal-noise-driven magnetization dynamics in a synthetic antiferromagnet*. arXiv:2001.02403. Available at: <https://doi.org/10.48550/arXiv.2001.02403>.
12. Menarini, M. and Lomakin, V. (2020) 'Thermal fluctuations in the Landau-Lifshitz-Bloch model', *Physical Review B*, 102(2), p. 024428. Available at: <https://doi.org/10.1103/PhysRevB.102.024428>.
13. Miceli, R. and McGuigan, M. (2019) *Thermo field dynamics on a quantum computer*. arXiv:1911.03335. Available at: <https://doi.org/10.48550/arXiv.1911.03335>.
14. Mondal, P., Suresh, A. and Nikolić, B.K. (2021) 'When can localized spins interacting with conduction electrons in ferro- or antiferromagnets be described classically via the Landau-Lifshitz equation: Transition from quantum many-body entangled to quantum-classical nonequilibrium states', *Physical Review B*, 104(21), p. 214401. Available at: <https://doi.org/10.1103/PhysRevB.104.214401>.
15. Mølmer, K., Castin, Y. and Dalibard, J. (1993) 'Monte Carlo wave-function method in quantum optics', *Journal of the Optical Society of America B*, 10(3), pp. 524–538. Available at: <https://www.phys.ens.psl.eu/~dalibard/publi3/osa_93.pdf>.
16. Nieves, P., Serantes, D., Atxitia, U. and Chubykalo-Fesenko, O. (2014) 'Quantum Landau-Lifshitz-Bloch (Quantum LLB) equation and its comparison with the classical case', *Physical Review B*, 90(10), p. 104428. Available at: <https://doi.org/10.1103/PhysRevB.90.104428>.
17. 'Quantum jump method' (2025) *Wikipedia*. Available at: <https://en.wikipedia.org/wiki/Quantum_jump_method>.
18. Rau, C., Jin, C. and Robert, M. (1988) 'Ferromagnetic order at Tb surfaces above the bulk Curie temperature', *Journal of Applied Physics*, 63(8), pp. 3667–3668. Available at: <https://doi.org/10.1063/1.340051>.
19. Schlimgen, A.W., Head-Marsden, K., Sager, L.M., Narang, P. and Mazziotti, D.A. (2202) 'Quantum simulation of the Lindblad equation using a unitary decomposition of operators', *Physical Review Research*, 4(2), p. 023216. Available at: <https://doi.org/10.1103/PhysRevResearch.4.023216>.
20. Sergi, A., Lamberto, D., Migliore, A. and Messina, A. (2023) 'Quantum–Classical Hybrid Systems and Ehrenfest's Theorem', *Entropy*, 25(4), p. 602. Available at: <https://doi.org/10.3390/e25040602>.
21. Su, V.P. (2021) 'Variational preparation of the thermofield double state of the Sachdev-Ye-Kitaev model', *Physical Review A*, 104(1), p. 012427. Available at: <https://link.aps.org/doi/10.1103/PhysRevA.104.012427>.
22. Thibaudeau, P., Fattouhi, M. and Buda-Prejbeanu, L.D. (2025) *Dynamic Landau-Lifshitz-Bloch-Slonczewski equations for spintronics*. [Preprint]. arXiv:2510.04562v2. Available at: <https://arxiv.org/abs/2510.04562v2>.
23. Wang, Y., Mulvihill, E., Hu, Z., Lyu, N., Shivpuje, S., Liu, Y., Soley, M.B., Geva, E., Batista, V.S. and Kais, S. (2022) *Simulating Open Quantum System Dynamics on NISQ Computers with Generalized Quantum Master Equations*. arXiv:2209.04956. Available at: <https://doi.org/10.48550/arXiv.2209.04956>.
24. Wieser, R. (2013) 'Comparison of Quantum and Classical Relaxation in Spin Dynamics', *Physical Review Letters*, 110(14), p. 147201. Available at: <https://doi.org/10.1103/PhysRevLett.110.147201>.
25. Wu, J. and Hsieh, T.H. (2019) 'Variational Thermal Quantum Simulation via Thermofield Double States', *Physical Review Letters*, 123(22), p. 220502. Available at: <https://link.aps.org/doi/10.1103/PhysRevLett.123.220502>.
