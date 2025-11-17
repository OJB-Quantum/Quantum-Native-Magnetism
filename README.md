# Quantum-Native-Magnetism
A repository for performing scientific quantum computing for atomistic and nanoscale magnetism studies. Created by Onri Jay Benally. 

[![License: GNU GPL v3](https://img.shields.io/badge/GNU_General_Public-License-Green)](https://choosealicense.com/licenses/gpl-3.0/)


---

Some useful up-front options are mentioned here and are later expanded upon in the this readme document. One can also think of this work as a way of applying the knowledge of various quantum-native solvers from the [Quantum-Native-Solvers](https://github.com/OJB-Quantum/Quantum-Native-Solvers) repository. 

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
    (accuracy referenced by a GPU Lindblad/kinetic integrator, scaled by a 156-qubit IBM Heron r2 processor with quantum error mitigation)

├─ 0. Classical baseline
│  └─ LLB–Boltzmann
│     ├─ LLB for macrospin m(t), transverse and longitudinal damping
│     ├─ Boltzmann for electrons and phonons, ultrafast kinetics and spin‑flip channels
│     └─ qLLB‑calibrated fallback:
│         replace α∥(T), α⊥(T), and m_eq(T) with qLLB‑consistent forms (KMS/detailed balance enforced)
│         whenever explicit quantum channels are unavailable, so the classical drift remains faithful near Tc
│
├─ 1. Quantum‑augmented effective models, physics backbone
│  ├─ qLLB (Nieves et al.) → quantum‑corrected rates and statistics, equilibrium m_eq^q(T)
│  ├─ Ehrenfest quantum drift:
│  │   • E‑qLLB‑B (steady regime), with Boltzmann kinetic envelope on the Bloch ball
│  │   • d‑E‑qLLB‑B (time‑dependent T(t), B_eff(t)), with quantum rates and equilibrium target
│  ├─ Quantum LLG (q‑LLG) and related forms, for cross‑checks and low‑damping limits
│  └─ Ehrenfest coupling to classical micromagnetic fields, plus kinetic collisions and diffusion
│
├─ 2. Dual‑track solver layer
│  ├─ 2A. GPU accuracy anchor, primary for Dy or Tb 1 nm Curie‑T extraction
│  │   ├─ Qiskit Dynamics LindbladModel, time‑dependent Hamiltonian and dissipators, with JAX/GPU acceleration
│  │   ├─ Implement qLLB‑consistent channels (KMS‑compliant T1/T2), finite‑volume Bloch‑ball Boltzmann
│  │   ├─ Compute m(T), χ(T), and Tc with finite‑size scaling, produce deterministic, mitigation‑free references
│  │   └─ Export MuMax3‑compatible kernels and lookup tables: α∥(T), α⊥(T), m_eq^q(T), diffusion tensors
│  │
│  └─ 2B. NISQ‑compatible quantum solvers (156‑qubit IBM Heron r2 + QEM via Qiskit Runtime)
│      ├─ Hardware + QEM envelope
│      │   ├─ Processor: heavy‑hex Heron r2, 156 qubits, tunable couplers, improved error rates
│      │   └─ QEM toolkit: M3 measurement mitigation, ZNE, Pauli twirling and randomized compiling,
│      │       selective PEC, hybrid mitigation passes integrated in Runtime
│      │
│      ├─ 2B‑1. QITE‑based E‑qLLB‑B/ d‑E‑qLLB‑B (primary QPU module)
│      │   ├─ Encode local E‑qLLB‑B generator on spin + bath registers, drawn from the 156‑qubit pool
│      │   ├─ Use VarQITE for imaginary‑time anchoring, add short real‑time segments for ⟨m(t)⟩ and ⟨E(t)⟩
│      │   ├─ Compile to heavy‑hex, optimize depth; apply ZNE + M3, add PEC on small subsystems when beneficial
│      │   └─ Validate against the GPU anchor; escalate only meshes or clusters that exceed GPU memory
│      │
│      ├─ 2B‑2. SSE‑based trajectories (quantum‑jump unraveling)
│      │   ├─ Variational SSE for local quantum spin clusters, tens of spins with ancillas
│      │   ├─ Optional GQME kernels for non‑Markovian baths and memory effects
│      │   └─ Apply QEM to ensemble expectation values; use as targeted, high‑fidelity subroutines
│      │
│      ├─ 2B‑3. QMS (Quantum‑Metropolis/ Metropolis‑Hastings) sampling
│      │   ├─ Prepare Gibbs states for Heisenberg + anisotropy models relevant to micromagnetics
│      │   ├─ Utilize the 156‑qubit space for moderate lattices that stress classical Monte Carlo
│      │   └─ Apply QEM to thermal observables, magnetization curves m(T) and susceptibilities χ(T)
│      │
│      └─ 2B‑4. TFD‑based finite‑T simulation
│          ├─ Variational TFD preparation on two copies, using up to ~1/2 of qubits per copy
│          ├─ Real‑time evolution of TFD states, finite‑temperature correlators and dynamical susceptibilities
│          └─ QEM during preparation and evolution, focusing on thermodynamic and linear‑response observables
│
└─ 3. Future fault‑tolerant and post‑NISQ directions
   ├─ Larger q‑LLG/q‑LLB clusters with exact SSE or full GKSL master‑equation solvers,
   │  across multiple Heron‑class nodes or on logical qubits with error correction
   ├─ Full QMS phase diagrams for realistic 2D/3D spin Hamiltonians, true quantum criticality and Curie‑T benchmarks
   └─ Hierarchical embedding
       ├─ Quantum clusters on logical qubits that learn qLLB/kinetic kernels
       └─ Classical MuMax3‑like grids that ingest those kernels at device scale
```

## Nanometer-scale Curie-T prediction

```
├─ Start: Ehrenfest-Bloch-ball models (MuMax3-friendly)
│  ├─ E-LLB-B  → Upgrade rates & m_eq → E-qLLB-B
│  └─ d-E-qLLB-B (time-dependent)  ← recommended practical target
│       └─ GPU Dynamics+JAX (primary), QPU QITE+QEM (scale)
└─ Conceptual high fidelity
   ├─ q-dLLB (GKSL) → equilibrium & dynamics
   └─ q-dLLB-B (GKSL + kinetic) → most complete
       └─ SSE/ QMS/ TFD as subroutines (GPU first, QPU selectively)
```

## Spin-thermal modeling

```
├─ Classical LLB (+Boltzmann) → E‑LLB‑B
│   └─ +quantum rates → E‑qLLB‑B
├─ Quantum LLB drift (+Boltzmann) → q‑LLB‑B 
│   ├─ Lower cost than full GKSL
│   ├─ KMS-detailed-balance rates, Brillouin m_eq^q(T)
│   └─ Great for Dy/Tb nanoscale Tc scans, moderate ramps
└─ Full GKSL density matrix
    ├─ q‑dLLB
    └─ +Boltzmann envelope → q‑dLLB‑B (gold standard, heaviest)
```

Note: Terms defined below.

---

For getting the Curie temperature of a nanometer-scale cuboid of relevant metal as accurately as possible, you should, very likely, run the quantum‑corrected models (qLLB‑based) on a GPU Lindblad solver because it integrates the open‑system equations directly and deterministically. Then, very usefully, you can prototype a MuMax3‑compatible “quantum module” on a Heron‑class QPU with QITE + QEM to reach bigger embedded problems, although it will be noisier than the GPU reference.

---

## A good start

| Technique | Drift structure  | Stationary thermodynamics  | Critical regime near ($T_c$) | Kinetic inputs  | Use case  |
|-|-|-|-|-|-|
| Ehrenfest–Landau‑Lifshitz–Boltzmann (E‑LL‑B)               | Precession + transverse relaxation ($\Gamma_\perp$) only; no longitudinal ( $\Gamma_\parallel$ ); magnitude conserved.                                                                                                                                     | Orientation‑only equilibrium on a fixed‑magnitude shell; no intrinsic route to ( $m \to 0$ ).                                                            | Misses magnitude collapse; cannot reproduce ($T_c$) without adding a longitudinal term; Binder cumulants ill‑posed.                                                                                            | Collision integrals may model directional diffusion and dephasing; still no critical behavior without ( $\Gamma_\parallel$ ).                                              | Directional diffusion, spin‑torque orientation control, domain‑wall transport far from ($T_c$); not reliable for ($T_c$).                                    |
| Ehrenfest–Landau‑Lifshitz‑Bloch–Boltzmann (E‑LLB‑B)             | Precession + ( $\Gamma_\perp$ ) + longitudinal ( $\Gamma_\parallel$ ) toward ( $m_{\mathrm{eq}}(T)$ ).                                                                                                                                                     | Gibbs/Boltzmann fixed point when rates satisfy KMS detailed balance; ( $m_{\mathrm{eq}}(T)!\to!0$ ) at ($T_c$).                                          | Captures magnitude collapse; supports ( $\chi(T)$ ) peaks and Binder cumulants; finite‑size scaling is possible.                                                                                               | Magnon/electron/phonon collision integrals; fluctuation–dissipation enforced; material‑specific rate models.                                                               | Curie‑temperature prediction, finite‑size scaling, and non‑equilibrium transients with classical relaxation statistics.                                      |
| Ehrenfest–quantum‑Landau‑Lifshitz‑Bloch–Boltzmann (E‑qLLB‑B)    | As E‑LLB‑B, but with quantum‑corrected ( $\Gamma_{\parallel,\perp}^{q}(T)$ ) and ( $m_{\mathrm{eq}}^{q}(T)$ ) from GKSL/KMS.                                                                                                                               | Quantum thermal fixed point; correct low‑($T$) saturation and near‑critical behavior under KMS.                                                          | Improved fidelity near ($T_c$), especially for 4f rare‑earths; reproduces slower quantum relaxation; robust cumulants.                                                                                         | Same structure, with quantum‑consistent fluctuation–dissipation; energy‑resolved scattering kernels if needed.                                                             | More accurate ($T_c$) for terbium at 1 nm; equilibrium scans and quasi‑static ramps with realistic rare‑earth statistics.                                    |
| Quantum Landau‑Lifshitz‑Bloch–Boltzmann (q‑LLB‑B)           | qLLB drift from quantum spin statistics: precession + ( $\Gamma_\perp^{q}$ ) + longitudinal ( $\Gamma_\parallel^{q}$ ) toward ( $m_{\mathrm{eq}}^{q}(T)$ ); magnitude not conserved; rates obey KMS; no full GKSL density matrix (that is q‑dLLB). | Quantum Gibbs fixed point with Brillouin ( $m_{\mathrm{eq}}^{q}(T)$ ); ( $m!\to!0$ ) at ($T_c$); correct low‑($T$) saturation; FDT‑consistent noise. | Captures magnitude collapse and quantum‑slowed relaxation; supports ( $\chi(T)$ ) peaks and Binder cumulants; finite‑size smearing via kinetic envelope; mean‑field‑like unless exchange/geometry refined. | Quantum‑consistent magnon/electron/phonon collision integrals; energy‑resolved kernels; KMS‑respecting stochastic driving; optional inter‑sublattice exchange pumping. | High‑fidelity ($T_c$) and ramps for nanoscale Tb/Dy when q‑dLLB‑B is too heavy; calibrates E‑qLLB‑B; equilibrium and moderately non‑adiabatic protocols. |
| Dynamic Ehrenfest–Landau‑Lifshitz‑Bloch–Boltzmann (d‑E‑LLB‑B)   | Time‑dependent ( $\mathbf{B}_ {eff}(t)$ ), ( $T(t)$ ); classical ( $\Gamma_{\parallel,\perp}(T(t))$ ) drive ( $m(t)$ ).                                                                                                                                    | Quasi‑static manifolds under slow ramps; no exact fixed point during drive; recovers E‑LLB‑B at steady state.                                            | Models dynamic crossing of ($T_c$) and ultrafast demagnetization; classical rates can bias Tb near‑critical behavior.                                                                                          | Time‑dependent collision integrals; non‑equilibrium FDT; optional spin‑transfer torque terms for driven devices.                                                           | Heat‑assisted writing, pump–probe demagnetization, and protocol design when quantum corrections are less dominant.                                           |
| Dynamic Ehrenfest–quantum-Landau‑Lifshitz‑Bloch–Boltzmann (d‑E‑qLLB‑B) | Time‑dependent quantum‑corrected ( $\Gamma_{\parallel,\perp}^{q}(T(t))$ ), ( $m_{\mathrm{eq}}^{q}(T(t))$ ), plus precession.                                                                                                                               | Tracks instantaneous quantum thermal targets under slow drives; supports controlled non‑adiabatic corrections.                                           | Best practical for dynamic ($T_c$) extraction with finite‑size and surface effects; accurate ramps and susceptibilities.                                                                                       | Time‑dependent, quantum‑consistent collision integrals; magnon–electron–phonon coupling; optional Slonczewski STT.                                                         | Primary choice for nanoscale Dy or Tb ($T_c$) under ramps, for ( $m(T)$ ), ( $\chi(T)$ ), and Binder analysis with kinetic broadening.                       |
| Quantum dynamic Landau‑Lifshitz‑Bloch (q‑dLLB)                  | Full GKSL master equation ( $\dot\rho=-\tfrac{i}{\hbar}[H,\rho]+\sum_\mu \mathcal D[L_\mu]\rho$ ); T1/T2 obey KMS; no Boltzmann layer.                                                                                                                          | Exact Gibbs fixed point for the specified model; ( $m!\to!0$ ) emerges when Hamiltonian + bath allow ordering.                                           | Captures amplitude collapse and critical response on small clusters; lacks kinetic‑ensemble broadening by itself.                                                                                              | Kinetics not explicit; approximate via disorder/size averaging, auxiliary reservoirs, or stochastic unraveling (SSE).                                                      | High‑fidelity cluster studies, equilibrium and short‑time dynamics; calibration of effective rates for Ehrenfest models.                                     |
| Quantum dynamic Landau‑Lifshitz‑Bloch–Boltzmann (q‑dLLB‑B)      | GKSL + kinetic Boltzmann envelope on the Bloch ball; most complete open‑system/kinetic coupling.                                                                                                                                                           | Quantum thermal steady state with kinetic consistency; detailed balance enforced both microscopically and macroscopically.                               | Captures magnitude collapse and critical fluctuations with finite‑size smearing; conceptually “gold standard.”                                                                                                 | Full, quantum‑consistent collision integrals; energy‑resolved, material‑specific kernels; strict fluctuation–dissipation.                                                  | Benchmarking and validation when computationally feasible; ultimate reference for nanoscale ($T_c$) and dynamic protocols.                                   |

---

## A more comprehensive approach for model selection

| Equation/ Approach | Quantum‑Native? (Q‑N) | What it captures near Tc | GPU path (Qiskit Dynamics/ Aer GPU) - Fit & Today’s Accuracy | QPU path (Heron 156q + QITE + QEM) - Fit & Today’s Accuracy | Best for 2 nm Dy or Tb Curie T? |
|-|-|-|-|-|-|
| E‑LLB‑B (Ehrenfest–Landau‑Lifshitz‑Bloch–Boltzmann)                     | No (classical); Q‑N‑compatible via GKSL emulation of T1/T2                    | Longitudinal + transverse damping, but with **classical** rates; kinetics via Boltzmann; no explicit quantum statistics                                        | **Good fit** via LindbladModel with effective rates; **Medium–Low accuracy** for Dy or Tb (misses quantum corrections)     | QITE+QEM can emulate effective channels, but depth/mitigation cost add bias; **Low–Medium accuracy**            | Baseline only                                               |
| E‑qLLB‑B (Ehrenfest–quantum‑Landau‑Lifshitz‑Bloch–Boltzmann)            | Q‑N‑compatible (rates & $m_eq$ from qLLB)                                       | Longitudinal collapse + **quantum‑corrected rates** and $m_eq^q(T)$; kinetics via Boltzmann; quasi‑static                                                        | **Very good fit** (native GKSL + JAX GPU); **High accuracy** for quasi‑static Tc extraction with finite‑size               | QITE(+TFD) + QEM workable, yet noise/mitigation bring variance; **Medium**                                      | Strong contender                                            |
| q‑LLB‑B (Quantum Landau‑Lifshitz‑Bloch–Boltzmann) | Q‑N‑compatible (qLLB drift + KMS‑balanced rates; not full density matrix) | **Quantum‑statistical** longitudinal & transverse channels with Brillouin ( $m_{\mathrm{eq}}^{q}(T)$ ); **kinetic envelope** adds finite‑size/transport broadening | **Very good–Excellent fit** (qLLB drift + GPU ODE/PDE; cheaper than q‑dLLB‑B); **High–Very High accuracy** for Dy/Tb ( $T_c$ ) | Use **small q‑dLLB blocks** to calibrate rates, integrate classically; **Medium–High** if calibration is stable | **Very strong** (just below q‑dLLB‑B for ultimate fidelity) |
| d‑E‑LLB‑B (Dynamic Ehrenfest–Landau‑Lifshitz‑Bloch–Boltzmann)           | No (classical); Q‑N‑compatible via GKSL                                       | Time‑dependent $T(t)$, $B_eff(t)$, classical rates; kinetics; handles ramps/pulses                                                                                 | **Good fit**, but **Medium** accuracy near Tc for Tb (classical bias)                                                      | Depth grows with steps; mitigation overhead; **Low–Medium**                                                     | Use for baselines                                           |
| d‑E‑qLLB‑B (Dynamic Ehrenfest–quantum‑Landau‑Lifshitz‑Bloch–Boltzmann)  | Q‑N‑compatible (qLLB‑derived)                                                 | **Quantum‑corrected** longitudinal & transverse channels **with time dependence**, plus kinetics; ideal for ramps                                              | **Best practical fit** (direct GKSL + JAX on GPU); **Very High accuracy** for Tc under ramps & finite‑size                 | Feasible with QITE(+TFD)/SSE + QEM; **Medium–High** if shallow; falls with step count                           | **Best (practical)**                                        |
| q‑dLLB (Quantum dynamic Landau‑Lifshitz‑Bloch, GKSL density‑matrix)     | **Yes (Q‑Native)**                                                            | Full GKSL Lindblad dynamics with **quantum statistics**, but **no kinetic f(m)** envelope                                                                      | **Excellent** for small/meso systems; **High** accuracy if finite‑size enters via Hamiltonian/parameters                   | **Best QPU candidate** (pure GKSL mapped to circuits), yet still mitigation‑limited; **Medium–High**            | Strong (needs size model)                                   |
| q‑dLLB‑B (Quantum dynamic Landau‑Lifshitz‑Bloch–Boltzmann)              | **Yes (Q‑Native)** (+ kinetic extension)                                      | As above **plus** kinetic distribution on Bloch ball; most complete near Tc, finite‑size, surfaces                                                             | **Gold standard** conceptually; heavy but **Very High** accuracy when implemented                                          | Very heavy (kinetic sampling + ancillas); **Medium** at best today                                              | **Best (if feasible)**                                      |

---

## Focused approach for a better "cost function" score

| Equation/ Approach  | Quantum‑Native? (Q‑N) | What it captures near Tc | GPU path (Qiskit Dynamics/ Aer GPU) - Fit & Today’s Accuracy| QPU path (Heron 156q + QITE + QEM) – Fit & Today’s Accuracy | Best for 2 nm Dy or Tb Curie T? |
|-|-|-|-|-|-|
| d‑E‑qLLB‑B (Dynamic Ehrenfest–quantum‑Landau‑Lifshitz‑Bloch–Boltzmann)  | Q‑N‑compatible (qLLB‑derived)               | **Quantum‑corrected** longitudinal & transverse channels **with time dependence**, plus kinetics; ideal for ramps   | **Best practical fit** (direct GKSL + JAX on GPU); **Very High accuracy** for Tc under ramps & finite‑size | Feasible with QITE(+TFD)/SSE + QEM; **Medium–High** if shallow; falls with step count                | **Best (practical)**            |
| q‑LLB‑B (Quantum Landau‑Lifshitz‑Bloch–Boltzmann) |Q‑N‑compatible (qLLB drift + KMS rates) | **Quantum‑slowed** relaxation + Brillouin ( $m_{\mathrm{eq}}^{q}(T)$ ); **kinetic** broadening; robust Binder cumulants | **High accuracy at moderate cost**; ideal when q‑dLLB‑B is too heavy                                       | **Use small GKSL sub‑blocks** to tune rates; main integration classical; **Medium–High**             | **Excellent trade‑off**         |
| q‑dLLB (Quantum dynamic Landau‑Lifshitz‑Bloch, GKSL density‑matrix)     | **Yes (Q‑Native)**                          | Full GKSL Lindblad dynamics with **quantum statistics**, but **no kinetic f(m)** envelope                           | **Excellent** for small/meso systems; **High** accuracy if finite‑size enters via Hamiltonian/parameters   | **Best QPU candidate** (pure GKSL mapped to circuits), yet still mitigation‑limited; **Medium–High** | Strong (needs size model)       |
| q‑dLLB‑B (Quantum dynamic Landau‑Lifshitz‑Bloch–Boltzmann)              | **Yes (Q‑Native)** (+ kinetic extension)    | As above **plus** kinetic distribution on Bloch ball; most complete near Tc, finite‑size, surfaces                  | **Gold standard** conceptually; heavy but **Very High** accuracy when implemented                          | Very heavy (kinetic sampling + ancillas); **Medium** at best today                                   | **Best (if feasible)**          |

## Tiny comparison table (what to run where, for nanometer ($T_C$))

| Task/ Form                                    | GPU (Qiskit Dynamics + JAX)                                                                                     | QPU (Heron-class + QITE + QEM)                                                                            |
|-|-|-|
| E‑LLB‑B (classical)                           | Fast baseline; not preferred near Tc (Tb)                                                                       | Feasible; use only as a scaffold; upgrade to qLLB asap                                                    |
| E‑qLLB‑B (quantum‑corrected)                  | Preferred baseline; steady‑state Tc scans                                                                       | Feasible; mitigated; validate vs GPU                                                                      |
| q‑LLB‑B (quantum‑native drift + kinetics)     | High accuracy at moderate cost; favored when q‑dLLB‑B is too heavy; calibrate vs small q‑dLLB or experiment     | Use QPU only to calibrate GKSL/KMS rates; main evolution classical; good anchor for mitigated studies     |
| d‑E‑qLLB‑B (time‑dependent)                   | Best practical for ramps + kinetics                                                                             | Feasible; deeper circuits & mitigation overhead                                                           |
| q‑dLLB                                        | High‑fidelity GKSL dynamics                                                                                     | Feasible on small subsystems; great as an anchor/check                                                    |
| q‑dLLB‑B                                      | Gold‑standard (GKSL + Boltzmann kinetics)                                                                       | Heavy today; reserve for focused, equilibrium subroutines                                                 |

**Why the GPU first:** Direct **Lindblad** integration with time‑dependent operators is native in Qiskit Dynamics and accelerates with JAX/GPU; mitigation‑free numerics dominate hardware noise for today’s long, dissipative evolutions. On the QPU, VarQITE/VarQRTE plus M3/ZNE/PEC is workable and valuable for scale/prototyping, but it carries residual bias/variance that grows with step count and circuit depth.

---

### **Near‑term approach:**

* **What to build:** Keep the MuMax3‑compatible module, but **promote the drift** from **E‑LLB‑B** to **E‑qLLB‑B** (and preferably **d‑E‑qLLB‑B** if you drive ( $T(t)$ ) or ( $\mathbf{B}* {\mathrm{eff}}(t))$ ). This swap keeps the Ehrenfest/Bloch‑ball structure, yet anchors the rates and ( $m* {\mathrm{eq}}(T)$ ) in a **quantum (GKSL/qLLB) derivation**, which is especially important near ( $T_C$ ) for rare‑earth terbium. On hardware you still run **QITE/VarQITE** with **QEM** (M3 measurement‑mitigation, ZNE extrapolation, and, when feasible, PEC). On day one you can stay within NISQ‑era depth limits on a **Heron‑class** system; in particular, the newer Heron‑family processors in IBM’s **System Two** environment provide improved error rates and scalable mitigation through the **Qiskit/Runtime** stack.

* **How to keep accuracy honest:** For the *same* model and discretization, run a **GPU reference** with **Qiskit Dynamics** (time‑dependent Lindblad) and **JAX/GPU** enabled. Use it to calibrate step sizes, validate QPU mitigation settings, and quantify dynamic‑ramp vs equilibrium bias when you extract ( $T_C$ ) from ( $m_{\mathrm{eq}}(T)$ ) or ( $\chi(T)$ ).

* **Why the upgrade from E‑LLB‑B?** Classical LLB rates mischaracterize quantum relaxation near ( $T_C$ ), whereas **qLLB** (quantum LLB) supplies **temperature‑dependent longitudinal/transverse rates** and the correct **( $m_{\mathrm{eq}}^{\mathrm{q}}(T)$ )**, which is exactly what matters in the narrow helical‑to‑ferromagnetic window of terbium (helical around ~231 K, ferromagnetic near ~219 K in bulk; nanoscale values are suppressed and broadened).

* **If something here “isn’t true,” how to *make it* true:** If you must start from **E‑LLB‑B** for code simplicity, **replace** the classical ( $\alpha_{\parallel,\perp}(T)$ ) and ( $m_{\mathrm{eq}}(T)$ ) by **qLLB‑consistent** rates and ( $m_{\mathrm{eq}}^{\mathrm{q}}(T)$ ) obtained from a GKSL model enforcing KMS detailed balance. At that point the model *becomes* **E‑qLLB‑B** (or **d‑E‑qLLB‑B**), i.e., quantum‑native compatible.

### **Conceptually best:**

* **What to favor conceptually:** **q‑dLLB** (full GKSL density‑matrix dynamics) and **q‑dLLB‑B** (adds the Boltzmann kinetic envelope) are the most faithful ways to model nanoscale Tb near ( $T_C$ ). For advanced observables—finite‑temperature entanglement, correlation functions, non‑Markovian baths—**SSE/quantum‑trajectories**, **QMS** (explicit master‑equations), and **TFD** (thermofield‑double thermal states) are the superior, general‑purpose tools.

* **How to use them now:** Implement them first on **GPU** (Dynamics for GKSL; Aer‑GPU or custom trajectories for SSE/TFD). Then, selectively port pieces to the **QPU** as **equilibrium‑only or small‑subsystem** subroutines with **VarQITE + QEM** (expect higher qubit counts, deeper circuits, larger shot budgets, and tougher optimization). In 2025 hardware terms, treat them as **high‑accuracy complements** to the Ehrenfest–qLLB–Boltzmann “main line” rather than replacements.


---

## Qiskit Aer GPU calculated results for Curie temperature prediction (among 4 candidate models)

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/5f76ed1b-d6f5-4e3c-8ec1-670107b9fb9f" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/58f36643-b78a-4848-9384-30543d5372f5" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/d54de531-13b5-4547-b996-41bd054ebfc7" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/73cc1dc2-2aca-4035-bb79-4776bdfeb4df" />

<img width="2177" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/9038e836-ac3e-4fc9-aaa3-f6082e7f4321" />

<img width="2177" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/3cf3f9d2-354a-43e0-9e84-9e1f752a5249" />

<img width="2177" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/49c19949-b76a-4f67-a18b-ed03a043e4db" />

<img width="2177" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/e952dfe6-444f-495b-b5b9-d6eb118425c8" />




## Qiskit Real QPU calculated results (with error mitigation) for Curie temperature prediction 

<img width="1433" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/7f75bea4-13e5-4c16-a7ed-d7e53d020655" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/dc9ce9da-9db3-4d52-a056-1ab03259449a" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/5bf449a7-6584-429c-8921-64abe70451c6" />

<img width="1386" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/6bdbe326-c117-4821-88e9-42dd8d6588bb" />


## Screenshot of the quantum circuit 

<img width="2461" height="125" alt="image" src="https://github.com/user-attachments/assets/088ae67e-515e-4261-b153-cd54df1502ce" />



---

## Relevant equations:

### Notation

* Magnetization (unit): ( $\boldsymbol{m}=\boldsymbol{M}/M_s$ ), ( $m=\lvert\boldsymbol{m}\rvert$ ).
* Effective field: ( $\boldsymbol{B}_{\mathrm{eff}}=-\partial \mathcal{F}/\partial \boldsymbol{m}$ ), with Zeeman, exchange, anisotropy, Dzyaloshinskii–Moriya (DMI), etc.
* Quantum expectation (site (i)): ( $\boldsymbol{s}_i(t)=\mathrm{Tr}!\big[\rho(t),\boldsymbol{\sigma}_i\big]$ ).
* Unit field direction: ( $\hat{\boldsymbol{b}}$ = $\boldsymbol{B}*{\mathrm{eff}}$ / $\lvert \boldsymbol{B} * {\mathrm{eff}} \rvert$ ).
* Rates: longitudinal ( $\Gamma_\parallel=1/T_1$ ), transverse ( $\Gamma_\perp=1/T_2$ ); LLB coefficients ( $\alpha_\parallel(T),\alpha_\perp(T)$ ).
* Equilibrium magnetization: ( $m_{\mathrm{eq}}(T)$ ) (or ( $m_{\mathrm{eq}}(T,D)$ ) if a finite-size parameter ( $D$ ) is included).
* Pauli vector: ( $\boldsymbol{\sigma}_i=(\sigma_i^x,\sigma_i^y,\sigma_i^z)$ ).
* Density matrix: ( $\rho(t)$ ); pure state trajectories: ( $\lvert\psi(t)\rangle$ ).
* Bohr frequencies: ( $\omega_{mn}=(E_m-E_n)/\hbar$ ), inverse temperature ( $\beta=1/k_{\mathrm B}T$ ).

---

## 0) Ehrenfest theorems (quantum → classical bridge)

### 0.1 Canonical Ehrenfest theorem (position–momentum pair)

For a particle of mass ( $m$ ) with Hamiltonian ( $H=\hat{p}^2/(2m)+V(\hat{x},t)$ ):

$$ \frac{d}{dt}\langle \hat{x}\rangle = \frac{1}{m},\langle \hat{p}\rangle, $$

$$ \frac{d}{dt}\langle \hat{p}\rangle = -\left\langle \frac{\partial V(\hat{x},t)}{\partial x}\right\rangle = \langle \hat{F}\rangle. $$

These are the canonical Ehrenfest relations: the expectation values obey Newton-like equations whenever the wave packet remains sufficiently localized.

### 0.2 General Ehrenfest theorem (any operator ( $\hat{A}$ ))

For a general (possibly time-dependent) operator ( $\hat{A}(t)$ ):

$$ \frac{d}{dt}\langle \hat{A}\rangle = \left\langle\frac{\partial \hat{A}}{\partial t}\right\rangle + \frac{1}{i\hbar},\big\langle[\hat{A},\hat{H}]\big\rangle. $$

This is the general Ehrenfest theorem, valid in Schrödinger or Heisenberg pictures. Choosing appropriate ( $\hat{A}$ ) (e.g., spin components, magnetization operators) yields all the Ehrenfest magnetization equations below.

---

## 1) Landau–Lifshitz / Gilbert and stochastic LLG

### 1.1 Classical original (Gilbert and LL forms)

**Gilbert form**:

$$ \frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}} + \frac{\alpha}{m},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt}. $$

**Landau–Lifshitz (LL) form**:

$$ \frac{d\boldsymbol{m}}{dt} = -\gamma',\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} - \gamma'\alpha,\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}}\big), \qquad \gamma' = \frac{\gamma}{1+\alpha^2}. $$

### 1.2 Stochastic LLG (thermal field)

$$ \frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\big(\boldsymbol{B}*{\mathrm{eff}}+\boldsymbol{B}*{\mathrm{th}}(t)\big) + \frac{\alpha}{m},\boldsymbol{m}\times\frac{d\boldsymbol{m}}{dt}. $$

$$ \big\langle B_{\mathrm{th},\mu}(t),B_{\mathrm{th},\nu}(t')\big\rangle = 2D,\delta_{\mu\nu},\delta(t-t'), \qquad D \propto \frac{\alpha k_{\mathrm B}T}{\gamma M_s V}. $$

This is the stochastic LLG underlying MuMax3-style micromagnetics.

### 1.3 Ehrenfest quantum Landau–Lifshitz (E-qLL)

Take a quantum spin (or spin cluster) with Hamiltonian ( $\hat{H}$ ) and define ( $\boldsymbol{s}(t)=\tfrac12\langle\boldsymbol{\sigma}(t)\rangle$ ). From the general Ehrenfest theorem:

$$ \frac{d\boldsymbol{s}}{dt} = \frac{1}{i\hbar},\big\langle[\tfrac12\boldsymbol{\sigma},\hat{H}]\big\rangle \equiv \boldsymbol{s}\times\boldsymbol{\Omega}_{\mathrm{eff}}. $$

Here ( $\boldsymbol{\Omega}_{\mathrm{eff}}$ ) is an operator-averaged precession vector (for ( $\hat{H}= -\tfrac{\hbar\gamma}{2}\boldsymbol{B}*{\mathrm{eff}}!\cdot!\boldsymbol{\sigma}$ one recovers the LL precession term).

Adding phenomenological transverse damping at the expectation level gives

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} - \Gamma*{\perp}\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big), \qquad \hat{\boldsymbol{b}} = \frac{\boldsymbol{B}*{\mathrm{eff}}}{\lvert\boldsymbol{B}*{\mathrm{eff}}\rvert}. $$

This is the Ehrenfest quantum Landau–Lifshitz (E-qLL) equation (transverse damping only).

### 1.4 Ehrenfest quantum Landau–Lifshitz–Gilbert (E-qLLG)

A quantum analog of LLG (q-LLG) can be formulated for the density operator ( $\rho$ ) of a spin system. At the level of ( $\boldsymbol{s}(t)=\mathrm{Tr}[\rho(t),\boldsymbol{\sigma}/2]$ ), one obtains an Ehrenfest-type qLLG:

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} + \alpha*{\mathrm{q}};\boldsymbol{s}\times\frac{d\boldsymbol{s}}{dt}. $$

Equivalently, in LL form,

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} - \frac{2\alpha*{\mathrm{q}}}{\hbar},\boldsymbol{s}\times\big(\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}\big). $$

Here ( $\alpha_{\mathrm{q}}$ ) is a quantum-derived damping coefficient constrained by the underlying q-LLG master equation. This is the Ehrenfest quantum Landau–Lifshitz–Gilbert (E-qLLG) equation.

### 1.5 Ehrenfest quantum Landau–Lifshitz–Gilbert–Slonczewski (E-qLLGS)

Include a Slonczewski spin-transfer torque (STT) term with polarization ( $\hat{\boldsymbol{p}}$ ):

$$ \left.\frac{d\boldsymbol{s}}{dt}\right|_{\mathrm{STT}} = -\gamma a_J,\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big). $$

The full Ehrenfest qLLG–Slonczewski (E-qLLGS) equation becomes

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} - \frac{2\alpha*{\mathrm{q}}}{\hbar},\boldsymbol{s}\times\big(\boldsymbol{s}\times\boldsymbol{B}_{\mathrm{eff}}\big) - \gamma a_J,\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big). $$

At the GKSL level this can be realized by adding spin-current-induced Lindblad channels or effective non-Hermitian terms coupled to a reservoir of polarized conduction electrons.

### 1.6 Adapted (Qiskit-ready, LL/LLG sector)

**Hamiltonian (single macrospin or per site (i))**

$$ \hat{H}(t) = -\frac{\hbar\gamma}{2},\boldsymbol{B}_{\mathrm{eff}}(t)\cdot\boldsymbol{\sigma} \quad \text{(encode as a sum of Pauli terms).} $$

**Optional LL/LLG-type damping on simulators**

Add appropriate Lindblad operators ( $L^\pm,L^z$ ) (see §4) so that the Ehrenfest equations reproduce E-qLL or E-qLLG. This is directly compatible with a Lindblad-based simulator or a QITE/ VQS approximation to Lindblad dynamics in Qiskit ≥ 2.2.

---

## 2) Landau–Lifshitz–Bloch (LLB) and dynamic extensions

### 2.1 Classical original LLB (single-sublattice macrospin)

$$ \frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\boldsymbol{B}* {\mathrm{eff}} + \frac{\gamma,\alpha*{\parallel}(T)}{m^2},\big(\boldsymbol{m}!\cdot!\boldsymbol{B}* {\mathrm{eff}}\big),\boldsymbol{m} - \frac{\gamma,\alpha*{\perp}(T)}{m^2},\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}\big). $$

A standard longitudinal field component inside ( $\boldsymbol{B}_{\mathrm{eff}}$ ) is

$$ \boldsymbol{B}* {\parallel} = \frac{1}{\chi*{\parallel}(T)} \left(1-\frac{m^2}{m_{\mathrm{eq}}^2(T)}\right)\boldsymbol{m}, \qquad \lvert \boldsymbol{m}\rvert \to m_{\mathrm{eq}}(T). $$

### 2.2 Dynamic LLB (time-dependent fields and parameters)

Allow explicit time dependence (e.g. ultrafast laser pumping, evolving temperature ( $T(t)$ )):

$$ \frac{d\boldsymbol{m}}{dt} = -\gamma,\boldsymbol{m}\times\boldsymbol{B}* {\mathrm{eff}}(\boldsymbol{m},t) + \frac{\gamma,\alpha*{\parallel}(T(t))}{m^2},\big(\boldsymbol{m}!\cdot!\boldsymbol{B}* {\mathrm{eff}}(\boldsymbol{m},t)\big),\boldsymbol{m} - \frac{\gamma,\alpha*{\perp}(T(t))}{m^2},\boldsymbol{m}\times\big(\boldsymbol{m}\times\boldsymbol{B}_{\mathrm{eff}}(\boldsymbol{m},t)\big). $$

This is the dynamic Landau–Lifshitz–Bloch (d-LLB) equation.

### 2.3 Dynamic LLB–Slonczewski (d-LLBS)

Add Slonczewski spin-transfer torque to d-LLB:

$$ \left.\frac{d\boldsymbol{m}}{dt}\right|_{\mathrm{STT}} = -\gamma a_J(t),\boldsymbol{m}\times\big(\boldsymbol{m}\times\hat{\boldsymbol{p}}\big). $$

So

$$ \frac{d\boldsymbol{m}}{dt} = \text{d-LLB drift from §2.2} - \gamma a_J(t),\boldsymbol{m}\times\big(\boldsymbol{m}\times\hat{\boldsymbol{p}}\big). $$

This is the dynamic Landau–Lifshitz–Bloch–Slonczewski (d-LLBS) equation.

### 2.4 Ehrenfest classical LLB (deterministic mean)

Replace ( $\boldsymbol{m}$ ) by quantum spin expectation ( $\boldsymbol{s}$ ) and write the LLB in Ehrenfest form:

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}*{\mathrm{eff}} - \Gamma* {\perp}(T)\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big) - \Gamma_{\parallel}(T)\Big[(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}}. $$

Here ( $\hat{\boldsymbol{b}}=\boldsymbol{B}* {\mathrm{eff}}/\lvert\boldsymbol{B}* {\mathrm{eff}}\rvert$ ) and classical LLB rates are related to ( $\alpha_{\parallel,\perp}(T)$ ).

### 2.5 Quantum Landau–Lifshitz–Bloch (qLLB) — density-matrix form

At the master-equation level:

$$ \dot{\rho} = -\frac{i}{\hbar}[\hat{H},\rho] + \sum_{\mu}\mathcal{D}[L_\mu]\rho, $$

with spin-flip and dephasing jump operators (e.g., ( $L^\pm\propto\sigma^\pm$ ), ( $L^z\propto\sigma^z$ )) chosen so that the Ehrenfest equation for ( $\boldsymbol{s}(t)=\mathrm{Tr}[\rho,\boldsymbol{\sigma}/2]$ ) reproduces an LLB-type structure, but with quantum-corrected rates:

$$ \alpha_{\parallel}(T) = \lambda,\frac{2T}{3T_C},\frac{2q_s}{\sinh(2q_s)}, \qquad \alpha_{\perp}(T) = \lambda\left[\frac{\tanh q_s}{q_s}-\frac{2T}{3T_C}\right], \qquad q_s = \frac{\mu H_{\mathrm{MFA}}}{k_{\mathrm B}T}. $$

### 2.6 Ehrenfest quantum LLB (E-qLLB)

From qLLB, the Ehrenfest equation for ( $\boldsymbol{s}$ ) yields Ehrenfest quantum LLB:

$$ \frac{d\boldsymbol{s}}{dt} = \frac{2}{\hbar},\boldsymbol{s}\times\boldsymbol{B}* {\mathrm{eff}} - \Gamma*{\perp}^{\mathrm{q}}(T)\Big(\boldsymbol{s}-(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big) - \Gamma_{\parallel}^{\mathrm{q}}(T)\Big[(\boldsymbol{s}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}^{\mathrm{q}}(T)\Big]\hat{\boldsymbol{b}}. $$

Here ( $\Gamma_{\parallel,\perp}^{\mathrm{q}}(T)$ ) and ( $m_{\mathrm{eq}}^{\mathrm{q}}(T)$ ) are obtained from qLLB.

### 2.7 Quantum dynamic LLB (q-dLLB) and E-q-dLLB

Allow explicit time dependence ( $T(t)$, ( $\boldsymbol{B}_{\mathrm{eff}}(t)$ ), etc.) in qLLB:

$$ \dot{\rho}(t) = -\frac{i}{\hbar}[\hat{H}(t),\rho(t)] + \sum_{\mu}\mathcal{D}[L_\mu(t)]\rho(t). $$

The Ehrenfest equation becomes

$$ \frac{d\boldsymbol{s}}{dt} = \text{E-qLLB drift with } T \to T(t), ; \boldsymbol{B}*{\mathrm{eff}} \to \boldsymbol{B}*{\mathrm{eff}}(t), $$

giving quantum dynamic LLB (q-dLLB) for ( $\rho$ ) and Ehrenfest quantum dynamic LLB (E-q-dLLB) for ( $\boldsymbol{s}$ ).

### 2.8 Quantum dynamic LLB–Slonczewski (q-dLLBS) and E-q-dLLBS

Add STT-type channels / Hamiltonian terms to q-dLLB:

$$ \dot{\rho}(t) = -\frac{i}{\hbar}\big[\hat{H}(t)+\hat{H}*{\mathrm{STT}}(t),\rho(t)\big] + \sum*{\mu}\mathcal{D}[L_\mu(t)]\rho(t), $$

with ( $\hat{H}_{\mathrm{STT}} \propto a_J(t),\hat{\boldsymbol{S}}\cdot(\hat{\boldsymbol{S}}\times\hat{\boldsymbol{p}})$ ) in an appropriate spin representation. At the Ehrenfest level this yields Ehrenfest quantum dynamic LLB–Slonczewski (E-q-dLLBS):

$$ \frac{d\boldsymbol{s}}{dt} = \text{E-q-dLLB drift of §2.7} - \gamma a_J(t),\boldsymbol{s}\times\big(\boldsymbol{s}\times\hat{\boldsymbol{p}}\big). $$

---

## 3) Kinetic-aware Ehrenfest–Boltzmann models on the Bloch ball

Introduce a distribution ( $f(\boldsymbol{m},t)$ ) (or ( $f(\boldsymbol{s},t)$ )) on the Bloch ball and write a drift–diffusion–collision kinetic equation:

$$ \partial_t f + \nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{A}(\boldsymbol{m},T,t),f\big] = \nabla_{\boldsymbol{m}}!\cdot!\big[\boldsymbol{D}(\boldsymbol{m},T,t),\nabla_{\boldsymbol{m}} f\big] + \mathcal{C}[f], $$

where ( $\mathcal{C}[f]$ ) encodes magnon/electron/phonon collisions and ( $\boldsymbol{D}$ ) satisfies fluctuation–dissipation.

### 3.1 Ehrenfest–LLB–Boltzmann (E-LLB-B)

Use LLB drift:

$$ \boldsymbol{A}* {\mathrm{LLB}}(\boldsymbol{m},T,t) = \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} - \Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big) - \Gamma_{\parallel}(T)\Big[(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}(T)\Big]\hat{\boldsymbol{b}}. $$

With ( $\Gamma_{\parallel,\perp}$ ) and ( $m_{\mathrm{eq}}(T)$ ) obtained from classical LLB, this defines the Ehrenfest–LLB–Boltzmann (E-LLB-B) equation.

### 3.2 Ehrenfest–quantum-LLB–Boltzmann (E-qLLB-B)

Replace classical LLB rates by qLLB rates (§2.5):

$$ \boldsymbol{A}* {\mathrm{qLLB}}(\boldsymbol{m},T,t) = \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}*{\mathrm{eff}} - \Gamma_{\perp}^{\mathrm{q}}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big) - \Gamma_{\parallel}^{\mathrm{q}}(T)\Big[(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}})-m_{\mathrm{eq}}^{\mathrm{q}}(T)\Big]\hat{\boldsymbol{b}}. $$

Then Ehrenfest quantum LLB–Boltzmann (E-qLLB-B) is the kinetic equation with ( $\boldsymbol{A}_{\mathrm{qLLB}}$ ).

### 3.3 Ehrenfest–LL–Boltzmann (E-LL-B) and E-qLL-B

Use LL drift only:

$$ \boldsymbol{A}* {\mathrm{LL}}(\boldsymbol{m},T,t) = \frac{2}{\hbar},\boldsymbol{m}\times\boldsymbol{B}* {\mathrm{eff}} - \Gamma_{\perp}(T)\Big(\boldsymbol{m}-(\boldsymbol{m}!\cdot!\hat{\boldsymbol{b}}),\hat{\boldsymbol{b}}\Big). $$

This lacks a longitudinal channel, so it cannot alone produce the amplitude collapse ( $m\to 0$ ) at ( $T_C$ ). With quantum-derived ( $\Gamma_{\perp}^{\mathrm{q}}(T)$ ) this becomes Ehrenfest quantum LL–Boltzmann (E-qLL-B).

### 3.4 Dynamic and Slonczewski variants (d-E-LLB-B, d-E-qLLB-B-S)

For the dynamic and Slonczewski cases:

* Dynamic Ehrenfest–LLB–Boltzmann (d-E-LLB-B): use ( $\boldsymbol{A}_{\mathrm{LLB}}(\boldsymbol{m},T(t),t)$ ).
* Dynamic Ehrenfest–quantum-LLB–Boltzmann (d-E-qLLB-B): use ( $\boldsymbol{A}_{\mathrm{qLLB}}(\boldsymbol{m},T(t),t)$ ).
* Dynamic Ehrenfest–quantum-LLB–Boltzmann–Slonczewski (d-E-qLLB-B-S): add STT drift ( $-\gamma a_J(t),\boldsymbol{m}\times(\boldsymbol{m}\times\hat{\boldsymbol{p}})$ ) to ( $\boldsymbol{A}_{\mathrm{qLLB}}$ ).

These are the Ehrenfest quantum Landau–Lifshitz–Bloch–Boltzmann, Ehrenfest quantum Landau–Lifshitz–Boltzmann, dynamic and Slonczewski-extended Boltzmann couplings.

---

## 4) Lindblad (GKSL) backbone and thermal anchoring (quantum-native)

### 4.1 GKSL master equation

For all quantum-native versions (qLL, qLLG, qLLB, q-dLLB, q-dLLBS), the density matrix obeys a GKSL equation:

$$ \dot{\rho} = -\frac{i}{\hbar}[\hat{H},\rho] + \sum_\mu\Big(L_\mu,\rho,L_\mu^\dagger - \tfrac12{L_\mu^\dagger L_\mu,\rho}\Big) \equiv -\frac{i}{\hbar}[\hat{H},\rho] + \sum_\mu \mathcal{D}[L_\mu]\rho. $$

With suitable ( $L_\mu$ ), Landau–Lifshitz damping can be derived directly from Lindbladian dissipation, establishing LL/LLG/LLB as Ehrenfest limits of GKSL dynamics.

### 4.2 Detailed balance (KMS) and ( $T_1/T_2$ ) dictionary

For a two-level splitting ( $\hbar\omega$ ) and spin-flip Lindblad operators:

$$ L^{-} = \sqrt{\gamma_{\downarrow}},\sigma^{-}, \qquad L^{+} = \sqrt{\gamma_{\uparrow}},\sigma^{+}, \qquad L^{z} = \sqrt{\gamma_{\phi}},\sigma^{z}. $$

Imposing KMS detailed balance:

$$ \frac{\gamma_{\uparrow}}{\gamma_{\downarrow}} = e^{-\beta\hbar\omega}, \qquad \Gamma_{\parallel} = \gamma_{\downarrow}+\gamma_{\uparrow}, \qquad \Gamma_{\perp} = \frac{\Gamma_{\parallel}}{2}+\gamma_{\phi}. $$

gives the thermal ( $T_1/T_2$ ) dictionary underlying LLB/qLLB.

### 4.3 Quantum-native mappings for QITE, SSE, and TFD

All the GKSL-level equations in §§1–3 can be mapped to quantum-circuit-compatible algorithms:

1. **QITE/ VQS for Lindblad**
   Treat GKSL as a first-order linear differential equation and apply Variational Quantum Simulation (VQS) or Quantum Imaginary-Time Evolution (QITE) extensions to open systems.
   *Implementation:*

   * Pauli-decompose ( $\hat{H}$ ) and ( $L_\mu$ ) (cf. §5).
   * Use Qiskit ≥ 2.2’s primitive-based VQS/QITE infrastructure with Quantum Error Mitigation (ZNE, MEM, etc.) on a 156-qubit Heron backend.

2. **Stochastic Schrödinger Equation (SSE) trajectories**
   Unravel GKSL into a stochastic Schrödinger equation:

$$ d\lvert\psi(t)\rangle = -\frac{i}{\hbar}\hat{H}* {\mathrm{eff}}(t)\lvert\psi(t)\rangle,dt + \sum* \mu \Big(\frac{L_\mu}{\sqrt{p_\mu(t)}}-\mathbb{I}\Big)\lvert\psi(t)\rangle,dN_\mu(t), $$

   where ( $\hat{H} * {\mathrm{eff}}$ ) includes non-Hermitian contributions, and ( $dN* \mu(t)$ ) are Poisson increments. Variational SSE implementations on NISQ devices have been developed for time-local master equations and non-Markovian baths.

3. **Thermofield-Double (TFD)-based VQAs**
   For finite-temperature qLLB/q-dLLB/q-dLLBS, one can represent thermal states as TFD states on a doubled Hilbert space and optimize a variational circuit to approximate ( $\lvert\mathrm{TFD}(\beta)\rangle$ ), then evolve unitarily.
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

$$ \hat{H} = -\frac{\hbar\gamma}{2}\sum_i \boldsymbol{B}* {\mathrm{eff},i}!\cdot!\boldsymbol{\sigma}*i - \sum* {\langle i,j\rangle}\big(J_x X_iX_j+J_y Y_iY_j+J_z Z_iZ_j\big) + \sum* {\langle i,j\rangle}\boldsymbol{D}_{ij}!\cdot!\big(\boldsymbol{\sigma}* i\times\boldsymbol{\sigma}* j\big) + \hat{H}* {\mathrm{ani}} + \hat{H}* {\mathrm{STT}}. $$

**Lindblad channels (per site (i))**

$$ L_i^{-} = \sqrt{\gamma_{\downarrow,i}},\sigma_i^{-}, \qquad L_i^{+} = \sqrt{\gamma_{\uparrow,i}},\sigma_i^{+}, \qquad L_i^{z} = \sqrt{\gamma_{\phi,i}},\sigma_i^{z}. $$

With KMS-consistent rates as in §4.2, these primitives reproduce Ehrenfest LL/LLG/LLB/qLLB drifts when inserted into the GKSL backbone (§4.1), and their Pauli decompositions are directly usable in Qiskit circuits.

---

## 6) Index of equations (chronological/ conceptual)

* Canonical Ehrenfest theorem: §0.1.

* General Ehrenfest theorem: §0.2.

* LLG/ LL (classical): §1.1.

* sLLG: thermal-field form + white-noise correlator, §1.2.

* Ehrenfest quantum LL (E-qLL): §1.3.

* Ehrenfest quantum LLG (E-qLLG): §1.4.

* Ehrenfest quantum LLG–Slonczewski (E-qLLGS): §1.5.

* Classical LLB: §2.1.

* Dynamic LLB (d-LLB): §2.2.

* Dynamic LLB–Slonczewski (d-LLBS): §2.3.

* Ehrenfest classical LLB: §2.4.

* Quantum LLB (qLLB, GKSL form): §2.5.

* Ehrenfest quantum LLB (E-qLLB): §2.6.

* Quantum dynamic LLB (q-dLLB) + E-q-dLLB: §2.7.

* Quantum dynamic LLB–Slonczewski (q-dLLBS) + E-q-dLLBS: §2.8.

* Ehrenfest–LLB–Boltzmann (E-LLB-B): §3.1.

* Ehrenfest–quantum-LLB–Boltzmann (E-qLLB-B): §3.2.

* Ehrenfest–LL–Boltzmann (E-LL-B) and E-qLL-B: §3.3.

* Dynamic Ehrenfest–(q)LLB–Boltzmann(+Slonczewski): §3.4.

* GKSL master equation (quantum-native backbone): §4.1.

* KMS + ( $T_1/T_2$ ) dictionary: §4.2.

* Quantum-native QITE/ SSE/ TFD mappings (Qiskit-compatible): §4.3.

* Adapted Hamiltonian + dissipators (Qiskit Pauli form): §5.

---

### Symbol mini‑glossary

( $\gamma$ ): gyromagnetic ratio; ( $\alpha$ ): Gilbert damping; ( $\alpha_{\parallel,\perp}(T)$ ): LLB coefficients; ( $\Gamma_{\parallel,\perp}$ ): Ehrenfest rates; ( $m_{\mathrm{eq}}(T)$ ): equilibrium magnetization; ( $\chi_{\parallel}(T)$ ): longitudinal susceptibility; ( $\boldsymbol{\sigma}$ = $(X,Y,Z)$ ): Pauli vector; ( $J_{x,y,z}$ ): exchange couplings; ( $\boldsymbol{D} * {ij}$ ): DMI vector; ( $\gamma*{\uparrow,\downarrow,\phi}$ ): Lindblad rate parameters; ( $\beta=(k_{\mathrm B}T)^{-1}$ ).
## Acronym glossary

```
LLB  : Landau–Lifshitz–Bloch (finite‑T magnetization dynamics)
qLLB : Quantum LLB (GKSL‑derived rates and m_eq^q(T))
E‑(q)LLB‑B : Ehrenfest (quantum) LLB + Boltzmann kinetics on Bloch ball
d‑E‑(q)LLB‑B : Dynamic E‑(q)LLB‑B with T(t), B_eff(t)
q‑dLLB : Quantum dynamic LLB (full GKSL density‑matrix)
q‑dLLB‑B : q‑dLLB with kinetic collisions/diffusion
GKSL : Gorini–Kossakowski–Sudarshan–Lindblad (quantum master equation)
QITE / VarQITE : (Variational) Quantum Imaginary‑Time Evolution
QEM : Quantum Error Mitigation (M3, ZNE, PEC)
M3 : Matrix‑free Measurement Mitigation
ZNE : Zero‑Noise Extrapolation
TFD : Thermofield Double (finite‑T pure‑state representation)
SSE : Stochastic Schrödinger Equation (quantum trajectories)
QMS : Quantum Master‑Equation solvers (e.g., GKSL integrators)
```


---

### **References**

1. Atxitia, U., Hinzke, D. and Nowak, U. (2017) 'Fundamentals and applications of the Landau–Lifshitz–Bloch equation', *Journal of Physics D: Applied Physics*, 50(3), p. 033003. Available at: <https://doi.org/10.1088/1361-6463/50/3/033003>.

2. Castin, Y., Dalibard, J. and Mølmer, K. (2008) *A Wave Function approach to dissipative processes*. arXiv:0805.4002. Available at: <https://doi.org/10.48550/arXiv.0805.4002>.

3. Chen, C-F., Kastoryano, M., Brandão, F.G.S.L. and Gilyén, A. (2025) 'Efficient quantum thermal simulation', *Nature*, 646(8085), pp. 561–566. Available at: <https://doi.org/10.1038/s41586-025-09583-x>.

4. Dalibard, J., Castin, Y. and Mølmer, K. (1992) 'Wave-function approach to dissipative processes in quantum optics', *Physical Review Letters*, 68(5), pp. 580–583. Available at: <https://doi.org/10.1103/PhysRevLett.68.580>.

5. Devyaterikov, D.I., Proglyado, V.V., Zhaketov, V.D., Nikitenko, Y.V., Kondrat'ev, O.A., Pashaev, E.M., Subbotin, I.A., Zverev, V.I., Kravtsov, E.A. and Ustinov, V.V. (2021) 'Influence of Dimensional Effects on the Curie Temperature of Dy and Ho Thin Films', *Physics of Metals and Metallography*, 122(5), pp. 465–471. Available at: <https://doi.org/10.1134/S0031918X21050033>.

6. Donahue, M.J. and Porter, D.G. (1999) *OOMMF User's Guide, Version 1.0*. Interagency Report NISTIR 6376. National Institute of Standards and Technology, Gaithersburg, MD. Available at: <https://doi.org/10.6028/NIST.IR.6376>.

7. Evans, R.F.L., Hinzke, D., Atxitia, U., Nowak, U., Chantrell, R.W. and Chubykalo-Fesenko, O. (2012) 'Stochastic form of the Landau-Lifshitz-Bloch equation', *Physical Review B*, 85(1), p. 014433. Available at: <https://doi.org/10.1103/PhysRevB.85.014433>.

8. Garanin, D.A. (1998) *Fokker-Planck and Landau-Lifshitz-Bloch equations for classical ferromagnets*. arXiv:cond-mat/9805054. Available at: <https://doi.org/10.48550/arXiv.cond-mat/9805054>.

9. Gokhale, S. and Manna, U. (2023) *Optimal control of the stochastic Landau-Lifshitz-Bloch equation*. arXiv:2305.10861. Available at: <https://doi.org/10.48550/arXiv.2305.10861>.

10. Gorini, V., Kossakowski, A. and Sudarshan, E.C.G. (1976) 'Completely positive dynamical semigroups of N-level systems', *Journal of Mathematical Physics*, 17(5), pp. 821–825. Available at: <https://doi.org/10.1063/1.522979>.

11. Liang, J-M., Lv, Q-Q., Wang, Z-X. and Fei, S-M. (2023) 'Assisted quantum simulation of open quantum systems', iScience, 26(4), p. 106306. Available at: <https://doi.org/10.1016/j.isci.2023.106306>.

12. Liu, S.H., Behrendt, D.R., Legvold, S. and Good, R.H., Jr. (1959) 'Interpretation of Magnetic Properties of Dysprosium', *Physical Review*, 116(6), pp. 1464–1468. Available at: <https://doi.org/10.1103/PhysRev.116.1464>.

13. Meil, D., Evelt, M., Pfau, B., Kläui, M., Atxitia, U. and Nowak, U. (2020) *Thermal-noise-driven magnetization dynamics in a synthetic antiferromagnet*. arXiv:2001.02403. Available at: <https://doi.org/10.48550/arXiv.2001.02403>.

14. Menarini, M. and Lomakin, V. (2020) 'Thermal fluctuations in the Landau-Lifshitz-Bloch model', *Physical Review B*, 102(2), p. 024428. Available at: <https://doi.org/10.1103/PhysRevB.102.024428>.

15. Miceli, R. and McGuigan, M. (2019) *Thermo field dynamics on a quantum computer*. arXiv:1911.03335. Available at: <https://doi.org/10.48550/arXiv.1911.03335>.

16. Mondal, P., Suresh, A. and Nikolić, B.K. (2021) 'When can localized spins interacting with conduction electrons in ferro- or antiferromagnets be described classically via the Landau-Lifshitz equation: Transition from quantum many-body entangled to quantum-classical nonequilibrium states', *Physical Review B*, 104(21), p. 214401. Available at: <https://doi.org/10.1103/PhysRevB.104.214401>.

17. Mølmer, K., Castin, Y. and Dalibard, J. (1993) 'Monte Carlo wave-function method in quantum optics', *Journal of the Optical Society of America B*, 10(3), pp. 524–538. Available at: <https://www.phys.ens.psl.eu/~dalibard/publi3/osa_93.pdf>.

18. Nieves, P., Serantes, D., Atxitia, U. and Chubykalo-Fesenko, O. (2014) 'Quantum Landau-Lifshitz-Bloch (Quantum LLB) equation and its comparison with the classical case', *Physical Review B*, 90(10), p. 104428. Available at: <https://doi.org/10.1103/PhysRevB.90.104428>.

19. 'Quantum jump method' (2025) *Wikipedia*. Available at: <https://en.wikipedia.org/wiki/Quantum_jump_method>.

20. Rau, C., Jin, C. and Robert, M. (1988) 'Ferromagnetic order at Tb surfaces above the bulk Curie temperature', *Journal of Applied Physics*, 63(8), pp. 3667–3668. Available at: <https://doi.org/10.1063/1.340051>.

21. Schlimgen, A.W., Head-Marsden, K., Sager, L.M., Narang, P. and Mazziotti, D.A. (2022) 'Quantum simulation of the Lindblad equation using a unitary decomposition of operators', *Physical Review Research*, 4(2), p. 023216. Available at: <https://doi.org/10.1103/PhysRevResearch.4.023216>.

22. Sergi, A., Lamberto, D., Migliore, A. and Messina, A. (2023) 'Quantum–Classical Hybrid Systems and Ehrenfest's Theorem', *Entropy*, 25(4), p. 602. Available at: <https://doi.org/10.3390/e25040602>.

23. Su, V.P. (2021) 'Variational preparation of the thermofield double state of the Sachdev-Ye-Kitaev model', *Physical Review A*, 104(1), p. 012427. Available at: <https://link.aps.org/doi/10.1103/PhysRevA.104.012427>.

24. Thibaudeau, P., Fattouhi, M. and Buda-Prejbeanu, L.D. (2025) *Dynamic Landau-Lifshitz-Bloch-Slonczewski equations for spintronics*. [Preprint]. arXiv:2510.04562v2. Available at: <https://arxiv.org/abs/2510.04562v2>.

25. Vansteenkiste, A., Leliaert, J., Dvornik, M., Hiasa, M., Garcia-Sanchez, F. and Van Waeyenberge, B. (2014) 'The design and verification of MuMax3', *AIP Advances*, 4(10), p. 107133. Available at: <https://doi.org/10.1063/1.4899186>.

26. Wang, Y., Mulvihill, E., Hu, Z., Lyu, N., Shivpuje, S., Liu, Y., Soley, M.B., Geva, E., Batista, V.S. and Kais, S. (2022) *Simulating Open Quantum System Dynamics on NISQ Computers with Generalized Quantum Master Equations*. arXiv:2209.04956. Available at: <https://doi.org/10.48550/arXiv.2209.04956>.

27. Wieser, R. (2013) 'Comparison of Quantum and Classical Relaxation in Spin Dynamics', *Physical Review Letters*, 110(14), p. 147201. Available at: <https://doi.org/10.1103/PhysRevLett.110.147201>.

28. Wu, J. and Hsieh, T.H. (2019) 'Variational Thermal Quantum Simulation via Thermofield Double States', *Physical Review Letters*, 123(22), p. 220502. Available at: <https://link.aps.org/doi/10.1103/PhysRevLett.123.220502>.
