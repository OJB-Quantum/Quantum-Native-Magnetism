# Quantum-Native-Magnetism
A repository for performing quantum-native scientific quantum computing for atomistic and nanoscale magnetism studies. Created by Onri Jay Benally. 

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
- QITE on a NISQ quantum computer implemented with quantum error mitigation on as many qubits as possible.
- Thermofield double-state on a NISQ or near-term quantum computer implemented with quantum error mitigation or quantum error correction.
- Stochastic Schrödinger‑equation-discretization-implemented near-term quantum computer with quantum error mitigation or quantum error correction.
- Quantum-metropolis-sampling-implemented near-term quantum computer with quantum error mitigation or quantum error correction.

---

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

## GPU calculated result for Curie temperature prediction 

<img width="1864" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/23c695e1-4b90-4a50-bbfd-f44ceb211da7" />


## Real QPU calculated result (without error mitigation) for Curie temperature prediction 

<img width="1864" height="auto" alt="Untitled" src="https://github.com/user-attachments/assets/b3f3bb24-bcda-4a7a-a094-70f9ddf6f22a" />


---

## Relevant equations: 

### Notation

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

with ($m_{\mathrm{eq}}$($T$)) from the Brillouin/Langevin mean‑field appropriate to the material.

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

---

References

Here is the reference list in numbered Harvard format, compatible with GitHub markdown.

The list has been alphabetized by author and duplicates from the original request have been consolidated. One entry (`https://arxiv.org/abs/2001.02403`) could not be verified and has been omitted.

1.  Atxitia, U., Hinzke, D. and Nowak, U. (2017) 'Fundamentals and applications of the Landau–Lifshitz–Bloch equation', *Journal of Physics D: Applied Physics*, 50(3), p. 033003. doi: 10.1088/1361-6463/50/3/033003.
2.  Castin, Y., Dalibard, J. and Mølmer, K. (2008) *A Wave Function approach to dissipative processes*. arXiv:0805.4002. doi: 10.48550/arXiv.0805.4002.
3.  Chen, C-F., Kastoryano, M., Brandão, F.G.S.L. and Gilyén, A. (2025) 'Efficient quantum thermal simulation', *Nature*, 646(8085), pp. 561–566. doi: 10.1038/s41586-025-09583-x.
4.  Dalibard, J., Castin, Y. and Mølmer, K. (1992) 'Wave-function approach to dissipative processes in quantum optics', *Physical Review Letters*, 68(5), pp. 580–583. doi: 10.1103/PhysRevLett.68.580.
5.  Devyaterikov, D.I., Proglyado, V.V., Zhaketov, V.D., Nikitenko, Y.V., Kondrat'ev, O.A., Pashaev, E.M., Subbotin, I.A., Zverev, V.I., Kravtsov, E.A. and Ustinov, V.V. (2021) 'Influence of Dimensional Effects on the Curie Temperature of Dy and Ho Thin Films', *Physics of Metals and Metallography*, 122(5), pp. 465–471. doi: 10.1134/S0031918X21050033.
6.  Evans, R.F.L., Hinzke, D., Atxitia, U., Nowak, U., Chantrell, R.W. and Chubykalo-Fesenko, O. (2012) 'Stochastic form of the Landau-Lifshitz-Bloch equation', *Physical Review B*, 85(1), p. 014433. doi: 10.1103/PhysRevB.85.014433.
7.  Garanin, D.A. (1998) *Fokker-Planck and Landau-Lifshitz-Bloch equations for classical ferromagnets*. arXiv:cond-mat/9805054. doi: 10.48550/arXiv.cond-mat/9805054.
8.  Gokhale, S. and Manna, U. (2023) *Optimal control of the stochastic Landau-Lifshitz-Bloch equation*. arXiv:2305.10861. doi: 10.48550/arXiv.2305.10861.
9.  Liu, S.H., Behrendt, D.R., Legvold, S. and Good, R.H., Jr. (1959) 'Interpretation of Magnetic Properties of Dysprosium', *Physical Review*, 116(6), pp. 1464–1468. doi: 10.1103/PhysRev.116.1464.
10. Meil, D., Evelt, M., Pfau, B., Kläui, M., Atxitia, U. and Nowak, U. (2020) Thermal-noise-driven magnetization dynamics in a synthetic antiferromagnet. arXiv:2001.02403. doi: 10.48550/arXiv.2001.02403. 
11. Menarini, M. and Lomakin, V. (2020) 'Thermal fluctuations in the Landau-Lifshitz-Bloch model', *Physical Review B*, 102(2), p. 024428. doi: 10.1103/PhysRevB.102.024428.
12. Miceli, R. and McGuigan, M. (2019) *Thermo field dynamics on a quantum computer*. arXiv:1911.03335. doi: 10.48550/arXiv.1911.03335.
13. Mondal, P., Suresh, A. and Nikolić, B.K. (2021) 'When can localized spins interacting with conduction electrons in ferro- or antiferromagnets be described classically via the Landau-Lifshitz equation: Transition from quantum many-body entangled to quantum-classical nonequilibrium states', *Physical Review B*, 104(21), p. 214401. doi: 10.1103/PhysRevB.104.214401.
14. Mølmer, K., Castin, Y. and Dalibard, J. (1993) 'Monte Carlo wave-function method in quantum optics', *Journal of the Optical Society of America B*, 10(3), pp. 524–538. Available at: [https://www.phys.ens.psl.eu/~dalibard/publi3/osa_93.pdf](https://www.phys.ens.psl.eu/~dalibard/publi3/osa_93.pdf).
15. Nieves, P., Serantes, D., Atxitia, U. and Chubykalo-Fesenko, O. (2014) 'Quantum Landau-Lifshitz-Bloch (Quantum LLB) equation and its comparison with the classical case', *Physical Review B*, 90(10), p. 104428. doi: 10.1103/PhysRevB.90.104428.
16. 'Quantum jump method' (2025) *Wikipedia*. Available at: [https://en.wikipedia.org/wiki/Quantum_jump_method](https://en.wikipedia.org/wiki/Quantum_jump_method) (Accessed: 12 November 2025).
17. Rau, C., Jin, C. and Robert, M. (1988) 'Ferromagnetic order at Tb surfaces above the bulk Curie temperature', *Journal of Applied Physics*, 63(8), pp. 3667–3668. doi: 10.1063/1.340051.
18. Schlimgen, A.W., Head-Marsden, K., Sager, L.M., Narang, P. and Mazziotti, D.A. (2022) 'Quantum simulation of the Lindblad equation using a unitary decomposition of operators', *Physical Review Research*, 4(2), p. 023216. doi: 10.1103/PhysRevResearch.4.023216.
19. Wieser, R. (2013) 'Comparison of Quantum and Classical Relaxation in Spin Dynamics', *Physical Review Letters*, 110(14), p. 147201. doi: 10.1103/PhysRevLett.110.147201.
