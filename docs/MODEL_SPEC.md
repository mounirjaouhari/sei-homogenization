# Model specification

Mathematical model implemented by the code in `src/`. Notation follows the paper.

## 1. Three-phase problem

A domain $\Omega\subset\mathbb R^d$ is split into electrode $\Omega_1^\varepsilon$,
SEI interphase $\mathcal B^\varepsilon$ (thickness $\varepsilon h$) and electrolyte
$\Omega_2^\varepsilon$. Quasi-static elasticity is coupled to transient lithium
diffusion through a compositional (Vegard) eigenstrain:

$$\nabla\!\cdot\boldsymbol\sigma^\varepsilon=\mathbf 0,\qquad
\boldsymbol\sigma^\varepsilon=\mathsf C^\varepsilon:\big(\boldsymbol\varepsilon(\mathbf u^\varepsilon)-\beta^\varepsilon(c^\varepsilon-c_0)\mathbf I\big),\qquad
\partial_t c^\varepsilon=\nabla\!\cdot(D^\varepsilon\nabla c^\varepsilon).$$

## 2. Scaling laws

$$\mathsf C^\varepsilon=\varepsilon^{\alpha}\widetilde{\mathsf C},\qquad
D^\varepsilon=\varepsilon^{\gamma}\widetilde D,\qquad
\beta^\varepsilon=\varepsilon^{\delta}\widetilde\beta,$$

with $\varepsilon$ the **dimensionless** thickness $\eta=\varepsilon/L$.

## 3. Effective jump conditions (matched asymptotics)

On the limit mid-surface $\mathcal S_0$, as $\varepsilon\to0$:

$$[\![\boldsymbol\sigma\!\cdot\!\mathbf n]\!]=\mathbf 0,\qquad
\boldsymbol\sigma\!\cdot\!\mathbf n=\mathsf K_{\mathrm{eff}}\!\cdot\!\big([\![\mathbf u]\!]-\boldsymbol\beta_{\mathrm{eff}}(\langle c\rangle-c_0)\big),\qquad
[\![c]\!]=R_{\mathrm{eff}}\langle J\rangle,\qquad[\![J]\!]=0.$$

## 4. Effective parameters (dimensionally-correct, FE-validated)

A layer of modulus $\mathsf C$, diffusivity $D$ and thickness $t=h\varepsilon$
behaves as a through-thickness spring and a diffusion resistance:

$$K_n=\frac{\tilde\lambda+2\tilde\mu}{t},\qquad
K_t=\frac{\tilde\mu}{t},\qquad
R_{\mathrm{eff}}=\frac{t}{D},\qquad
\beta_{\mathrm{eff}}=\beta\,t.$$

**Anisotropy** (validated in 2-D, untestable in 1-D):

$$\frac{K_n}{K_t}=\frac{\tilde\lambda+2\tilde\mu}{\tilde\mu}=\frac{2(1-\nu)}{1-2\nu}.$$

## 5. Asymptotic regimes (by $\alpha$)

| Regime | $\alpha$ | $\mathsf K_{\mathrm{eff}}$ | meaning |
|---|---|---|---|
| soft | $>1$ | $\to0$ | traction-free |
| intermediate | $=1$ | $O(1)$ | linear spring |
| critical | $=0$ | $\to\infty$, energy finite | cohesive law |
| stiff | $<0$ | $\to\infty$, energy $\to\infty$ | perfect bonding |

Diffusion regimes by $\gamma$: fast ($>1$), resistive ($=1$), blocking ($<1$).

## 6. Critical-regime cohesive law

$$\mathbf T=(1-d^\star)^2\,\widetilde{\mathsf K}\cdot\boldsymbol\Delta_{\mathrm{eff}},\qquad
\boldsymbol\Delta_{\mathrm{eff}}=[\![\mathbf u]\!]-\widetilde\beta\,h(\bar c-c_0)\widetilde{\mathsf A}\!\cdot\!(\widetilde{\mathsf C}\!:\!\mathbf I),$$

$$G_c=\varepsilon_{\mathrm{SEI}}\sum_i f_i g_{v,i},\qquad
\delta_c=\sqrt{2G_c/K_n}.$$

The volumetric fracture-energy density $g_{v,i}=\kappa\,\sigma_{r,i}^2/(2E_i)$ keeps
$G_c$ dimensionally consistent (units J/m²).

## 7. Convergence

The homogenized model converges to the resolved one at rate $\mathcal O(\varepsilon)$
(Theorem). Measured: super-linear $\sim\mathcal O(\varepsilon^{1.7\text{–}1.9})$ (1-D)
and $\mathcal O(\varepsilon^{0.97})$ for both normal and shear (2-D), for all
regimes.
