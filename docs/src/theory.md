# Theory

```@meta
CurrentModule = SaddlePaths
```

## Saddle Path Problems in Economics

### Mathematical Foundation

A saddle path problem arises in dynamic optimization when solving infinite-horizon optimal control problems. The canonical form is:

```math
\max_{c(t)} \int_0^{\infty} e^{-\rho t} u(c(t)) dt
```

subject to:
```math
\dot{k}(t) = f(k(t)) - c(t)
```
```math
k(0) = k_0 \text{ given}
```
```math
\lim_{t \to \infty} e^{-\rho t} \lambda(t) k(t) = 0 \text{ (transversality)}
```

### First-Order Conditions

The Hamiltonian for this problem is:
```math
H = u(c) + \lambda(f(k) - c)
```

The first-order conditions yield:
1. **Static efficiency**: $\frac{\partial H}{\partial c} = 0 \Rightarrow u'(c) = \lambda$
2. **Dynamic equation for costate**: $\dot{\lambda} = \rho \lambda - \frac{\partial H}{\partial k}$
3. **State equation**: $\dot{k} = f(k) - c$

### The Saddle Path

The system of differential equations:
```math
\begin{pmatrix}
\dot{k} \\
\dot{c}
\end{pmatrix} = 
\begin{pmatrix}
f(k) - c \\
c \cdot \frac{u'(c)}{u''(c)} \left( f'(k) - \rho \right)
\end{pmatrix}
```

has a steady state $(k^*, c^*)$ that is typically a saddle point:
- One stable eigenvalue (negative real part)
- One unstable eigenvalue (positive real part)

The **saddle path** is the unique trajectory converging to the steady state, corresponding to the stable manifold of the dynamical system.

## Classification of Variables

### State Variables

State variables represent stocks that:
- Accumulate over time
- Have initial conditions
- Cannot jump discontinuously

Examples:
- Physical capital stock ($k$)
- Human capital ($h$)
- Environmental quality ($e$)
- Debt levels ($b$)

### Costate Variables (Shadow Prices)

Costate variables represent:
- Shadow prices of state variables
- Control variables in the primal problem
- Variables that can jump at $t=0$

Examples:
- Consumption ($c$) as control
- Shadow price of capital ($\lambda$)
- Investment rate ($i$)
- Effort level ($e$)

### The Functional Notation

Writing $c(k)$ instead of $c$ indicates that consumption is a control variable that depends on the state. This notation:
1. Emphasizes the feedback nature of optimal control
2. Distinguishes controls from states
3. Aligns with the mathematical convention $c^*(k)$ for optimal policy

## Phase Diagram Analysis

### Two-Dimensional Systems

For a system with one state ($k$) and one costate ($c$):

```
    c
    ↑
    │     ╱ Saddle Path
    │    ╱
    │   ╱ 
    │  ● (k*, c*)
    │ ╱
    │╱____________→ k
```

The saddle path divides the phase space:
- Above: trajectories diverge (over-consumption)
- Below: trajectories violate transversality (under-consumption)
- On path: unique optimal trajectory

### Stability Analysis

Near the steady state, linearization gives:
```math
\begin{pmatrix}
\dot{k} \\
\dot{c}
\end{pmatrix} = 
J
\begin{pmatrix}
k - k^* \\
c - c^*
\end{pmatrix}
```

where the Jacobian $J$ has:
- Determinant: $\det(J) < 0$ (saddle point)
- Trace: $\text{tr}(J) = \lambda_1 + \lambda_2$ (sign ambiguous)

## Solution Methods

### 1. Shooting Method

1. Guess initial costate $c(0)$
2. Integrate forward the ODEs
3. Check if trajectory converges to steady state
4. Adjust $c(0)$ using bisection or Newton's method

### 2. Time Elimination

1. Solve for steady state
2. Use stable manifold theorem
3. Parameterize saddle path as $c = \psi(k)$
4. Solve the functional equation

### 3. Linearization Method

1. Linearize around steady state
2. Compute eigenvalues and eigenvectors
3. Stable eigenvector gives local saddle path
4. Use as initial approximation for global path

### 4. Projection Methods

1. Approximate policy function $c(k) = \sum_i a_i \phi_i(k)$
2. Use basis functions $\phi_i$ (Chebyshev, splines)
3. Solve for coefficients $a_i$ via collocation or Galerkin

## Common Economic Applications

### Growth Theory

**Ramsey-Cass-Koopmans Model**:
- State: capital stock
- Costate: consumption
- Saddle path: optimal consumption path

**Endogenous Growth**:
- States: physical capital, human capital, knowledge stock
- Costates: consumption, education time, R&D investment

### Environmental Economics

**Pollution Control**:
- States: capital, pollution stock
- Costates: consumption, abatement effort
- Trade-off: growth vs. environmental quality

### Resource Economics

**Hotelling Model**:
- State: resource stock
- Costate: extraction rate
- Saddle path: optimal depletion path

### Public Economics

**Optimal Taxation**:
- States: capital, government debt
- Costates: consumption, tax rates
- Saddle path: tax smoothing trajectory

## Numerical Considerations

### Stiff Systems

Economic models often yield stiff ODEs due to:
- Multiple time scales (fast consumption adjustment, slow capital accumulation)
- Wide parameter ranges
- Nonlinearities

Solutions:
- Implicit integration methods
- Adaptive time stepping
- Proper scaling of variables

### Boundary Conditions

The infinite horizon creates challenges:
- Initial condition: $k(0) = k_0$
- Terminal condition: transversality
- Saddle path: selecting $c(0)$ to satisfy both

### Convergence Criteria

For numerical solutions:
1. **Steady state convergence**: $|k(T) - k^*| < \epsilon$
2. **Transversality check**: $e^{-\rho T} \lambda(T) k(T) < \epsilon$
3. **Path stability**: eigenvalue analysis of linearized system

## Extensions

### Stochastic Saddle Paths

With uncertainty, the problem becomes:
```math
\max E_0 \left[ \int_0^{\infty} e^{-\rho t} u(c(t)) dt \right]
```

The saddle path becomes a stochastic process satisfying:
- Forward-looking: expectations of future states
- Backward-looking: current state realization

### Multiple Equilibria

Some models exhibit:
- Multiple steady states
- Multiple saddle paths
- Threshold effects and regime switching

### Time-Inconsistent Problems

With hyperbolic discounting or changing preferences:
- Saddle path depends on commitment ability
- Distinction between sophisticated and naive agents
- Game-theoretic interpretation (self vs. future selves)