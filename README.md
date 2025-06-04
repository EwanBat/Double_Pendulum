# Double Pendulum Simulation with 3D Stereographic Projection

## Author
**Ewan Bataille**

## Description

This project simulates the motion of a **double pendulum** using **stereographic projection** to represent the configuration on a plane. The simulation uses **Hamiltonian mechanics** to model the dynamics and numerically integrates the equations of motion to visualize the chaotic behavior of the system.

---

## 📚 Sections Overview

### 1. Stereographic Projection

We project the unit sphere (radius = 1, centered at \( O = (0,0,0) \)) onto the plane \( P: z = -1 \) from the north pole \( N = (0,0,1) \).

Any point $M = (\xi, \eta, \zeta)$ on the sphere satisfies: $\xi^2 + \eta^2 + \zeta^2 = 1 \tag{1.1}$

The stereographic projection from $N$ onto the plane gives coordinates: $x = \frac{2\xi}{1 - \zeta}, \quad y = \frac{2\eta}{1 - \zeta} \tag{1.2}$

The inverse transformation is: $(\xi, \eta, \zeta) = \left( \frac{2x}{1 + x^2 + y^2}, \frac{2y}{1 + x^2 + y^2}, \frac{-1 + x^2 + y^2}{1 + x^2 + y^2} \right) \tag{1.3}$

> ⚠️ The projection diverges at $$\( \zeta = 1 \)$$, so points near the north pole must be avoided in simulations.

---

### 2. Simple Pendulum Simulation

Using spherical coordinates, the motion of a pendulum of mass $\( m \)$ and length $\( l \)$ is governed by:

$$\
\begin{cases}
l \cdot \ddot{\theta} - l \cdot \dot{\varphi}^2 = -g \cdot \sin(\theta) \\
l \cdot \ddot{\varphi} + 2l \cdot \dot{\theta} \cdot \dot{\varphi} \cdot \cos(\theta) = 0
\end{cases}
\$$

This system can be solved numerically to obtain the 3D trajectory of the pendulum.

---

### 3. Double Pendulum with Hamiltonian Mechanics

We consider two pendulums with respective lengths $\( l_1, l_2 \)$ and masses $\( m_1, m_2 \)$. The first pendulum is anchored at the origin, and the second is attached to the end of the first.

Stereographic coordinates and momenta:

- For pendulum 1: $$\( q_1, q_2 \), \( p_1, p_2 \)$$
- For pendulum 2: $$\( q_3, q_4 \), \( p_3, p_4 \)$$

#### Hamiltonian Function

$$ H = \frac{h}{f} + l$$

$h = (A - B + C - D)^2 \cdot \frac{m_1}{m_2} + (A - B)^2 + \frac{m_2}{m_1} (y_1^2 + y_2^2) + (y_1 - y_3)^2 + (y_2 - y_4)^2 + 2F(q_1 y_3 + q_2 y_4) + 2G(q_3 y_1 + q_4 y_2) - 2FG(1 + q_1 q_3 + q_2 q_4)$

$f = 2 l_1^2 l_2^2 (m_2 a^2 b^2 + 4m_1 k(a b - k))$

$l = g((m_1 + m_2)(1 - \frac{2}{b})l_2 + m_1(1 - \frac{2}{a})l_1) \tag{2.1}$

- The first term represents **kinetic energy**
- The added term is **potential energy**

#### Definitions (from Appendix)

$g = 9.81 m/s²$ (gravitational constant)

$a = 1 + q1² + q2²$; $b = 1 + q3² + q4²$; $k = (q1 - q3)² + (q2 - q4)²$

$A = p1 * l2 * a * (a(q4 - q2) + q2 * k)$; $B = p2 * l2 * a * (a(q3 - q1) + q1 * k)$; $C = p3 * l1 * b * (b(q2 - q4) + q4 * k)$; $D = p4 * l1 * b * (b(q1 - q3) + q3 * k)$

$y1 = p1 * l2 * b * (b² / 2)$; $y2 = p2 * l2 * b * (b² / 2)$; $y3 = p3 * l1 * a * (b² / 2)$; $y4 = p4 * l1 * a * (b² / 2)$

$F = a * b * l2 * (q1 * p1 + q2 * p2)$; $G = a * b * l1 * (q3 * p3 + q4 * p4)$
