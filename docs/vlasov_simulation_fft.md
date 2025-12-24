# Solving the Vlasov-Poisson equation using Python

A plasma is a quasi-neutral gas of charged and neutral particles showing non-trivial collective behavior, distinct from other states of matter. The dynamics of plasma can be derived from the Boltzmann equation, where we can usually assume that we have a negligible amount of collisions. This leads us to the collisionless Boltzmann or **Vlasov equation**. In this tutorial, we want to simulate a simple one-dimensional plasma using the Vlasov equation and look at Landau damping and two kinetic instabilities. Since in one dimension we can't have any magnetic fields, we actually solve the electrostatic Vlasov-Poisson equation.

## Theory

If we want to simulate moving electrons, in principle we have to solve the full Maxwell system of equations, containing $E$ and $B$ fields. Since in the following we can safely assume the speed of the electrons $v$ to be much smaller than the speed of light $c$, i.e., $v\ll c$, we can neglect the contribution of the magnetic field and assume the effect of the electric field to be instantaneous. We know that the divergence of the electric field equals the charge density $$\nabla\cdot E=4\pi\rho$$ (in CGS units) and that in electrostatics the electric field can be derived from an electric potential $$E=-\nabla\phi.$$ Combined, we find Poissons equation for electrostatics $$\Delta\phi=4\pi\rho.$$ 

When simulating plasma, computational physicists usually use either the fluid discription of plasmas, called Magneto-Hydrodynmaics (MHD), or the semi-lagrangian Particle-In-Cell (PIC) method to sample the Vlasov equation. However, both approaches have their drawbacks. For example, the MHD theory is based on a few but fundamental assumptions which break down in certain, very interesting, regimes. PIC on the other hand does not have the limits on, e.g., resitivity, but needs a computationally very intensive amount of particles to generate solutions with a reasonable low amount of noise at a given space-time resolution.

The idea of simulating the Vlasov equation is to directly evolve the probability distribution function $f$ to find a charged particle in a given phase space region $\mathrm{d}x\mathrm{d}v$ at time $t$, i.e., the Vlasov equation is a 6+1 dimensional scalar function of position and velocity. However, since this is difficult to solve and we already learn a lot from simpler systems, we restrict the following discussion to the one-dimensional case. 

For a plasma there are usually two distribution functions, one for ions and one for electrons. But since ion mass is typically much larger than the electron mass ($m_i\gg m_e$), electrons more much faster than ions. In such a case, we can assume a neutralizing ion background and only simulate electron dynamics. Using our assumptions, the Vlasov equation can be written as $$\frac{\mathrm{d}}{\mathrm{d}t} f(x, v,t) = 0,$$
where charge density is then simply the first moment of the Vlasov equation times the electron charge $q_e$ $$\rho(x,t)=q_e\int f(x,v,t) \mathrm{d}v$$

Taking the total time derivative of $f$ using the chain rule we find the plasma kinetic equation $$\frac{\mathrm{d}}{\mathrm{d}t} f=\frac{\partial f}{\partial t}+\frac{\partial x}{\partial t}\frac{\partial f}{\partial x}+\frac{\partial v}{\partial t}\frac{\partial f}{\partial v}=0$$

By identifying $\frac{\partial x}{\partial t}=v$ and $\frac{\partial v}{\partial t}=a$ and using the equation for electrostatic force $F_E=m_ea=q_e E$ we derived the 1D Vlasov-Poisson equation: $$\frac{\partial f}{\partial t}+v\frac{\partial f}{\partial x}+\frac{q_\alpha E}{m}\frac{\partial f}{\partial v}=0.$$ The form of this equation resembles a _hyperbolic conservation law_ that can be computed numerically by employing standard techniques of fluid dynamics. For simplicity, in this tutorial we use the so-called upwind scheme.



# Implementation

To run the simulation of the Vlasov equation on a computer we have to implement a finite difference scheme. As our programming language of choice is python, we have to take special care on performance. This obscures the readability of the code but makes python code run nearly as fast as C code. First we need to import our favorite libraries. Here we will mainly rely on `numpy` for computation and `matplotlib` for visualization.



```python
import matplotlib.pyplot as plt
import numpy as np
import scipy
import tqdm
from maxpy.style import dark_mode

dark_mode(background="#191c1c")

CMAP = plt.get_cmap("cmr.chroma")
```

## Computing the electric field from a charge density using Fourier transforms

### Fourier transform conventions

We define the Fourier transform and its inverse as

$$
f(\mathbf{k}) = \int_{\mathbb{R}^3} f(\mathbf{r}) e^{-i\mathbf{k}\cdot\mathbf{r}} \, d^3 r
$$

$$
f(\mathbf{r}) = \frac{1}{(2\pi)^3} \int_{\mathbb{R}^3} f(\mathbf{k}) e^{i\mathbf{k}\cdot\mathbf{r}} \, d^3 k
$$

Different conventions may redistribute factors of $2\pi$ without changing physical results.

### Poisson’s equation in Fourier space

Using the correspondence

$$\nabla^2 \longrightarrow -k^2$$

Poisson’s equation becomes

$$-k^2 \phi(\mathbf{k}) = -4\pi\rho(\mathbf{k})$$

Solving for the potential gives

$$\phi(\mathbf{k}) = \frac{4\pi\rho(\mathbf{k})}{k^2}$$

This expression is valid for $k \neq 0$. The $k = 0$ mode corresponds to the total charge and requires separate treatment.

### Electric field in Fourier space

Since $\mathbf{E} = -\nabla \phi$ and

$$\nabla \longrightarrow i\mathbf{k}$$

the electric field in Fourier space is

$$\mathbf{E}(\mathbf{k}) = -i\mathbf{k}\,\phi(\mathbf{k})
= -i\mathbf{k}\,\frac{4\pi\rho(\mathbf{k})}{k^2}$$

or equivalently

$$\mathbf{E}(\mathbf{k}) = i\frac{\mathbf{k}}{k^2}4\pi\rho(\mathbf{k})$$

### Transforming back to real space

The real-space electric field is obtained by inverse Fourier transform:

$$\mathbf{E}(\mathbf{r}) =
\frac{1}{(2\pi)^3}
\int d^3 k \;
\frac{i\mathbf{k}}{k^2}\,
\rho(\mathbf{k})\,
e^{i\mathbf{k}\cdot\mathbf{r}}$$

This expression is equivalent to Coulomb’s law:

$$\mathbf{E}(\mathbf{r}) = \int d^3 r' \;
\rho(\mathbf{r}')
\frac{\mathbf{r}-\mathbf{r}'}{|\mathbf{r}-\mathbf{r}'|^3}$$

### Practical algorithm

Given $\rho(\mathbf{r})$:

1. Compute $\rho(\mathbf{k})$ via a Fourier transform.
2. Multiply by $\frac{i\mathbf{k}}{k^2}$ in Fourier space.
3. Apply the inverse Fourier transform to obtain $\mathbf{E}(\mathbf{r})$.

Using spectral solvers for the Poisson equation is especially useful when we deal with periodic boundary conditions, since they are satisfied by definition. The real benefit, however, comes from the Fast Fourier Transform (FFT).


```python
def solve_poisson_fft(rho, dx, nx):
    rho_k = np.fft.fft(rho)
    k_vals = 2 * np.pi * np.fft.fftfreq(nx, d=dx)

    E_k = np.zeros_like(rho_k, dtype=complex)
    nonzero = k_vals != 0
    E_k[nonzero] = (1j / k_vals[nonzero]) * rho_k[nonzero]

    E = np.real(np.fft.ifft(E_k))
    return E
```

## Upwind scheme

The upwind scheme is a finite difference or finite volume discretization method for hyperbolic partial differential equations. It accounts for the direction of information propagation.

### 1. Model problem

Consider the linear advection equation

$$\frac{\partial u}{\partial t} + a \frac{\partial u}{\partial x} = 0$$

where $a$ is a constant advection velocity. Characteristics propagate with speed $a$. The numerical flux is constructed using information from the upstream direction.

- If $a > 0$, information travels from left to right.
- If $a < 0$, information travels from right to left.

### First order upwind discretization

Let $u_i^n \approx u(x_i, t^n)$ with spatial step $\Delta x$ and time step $\Delta t$.

#### Case $a > 0$

$$\frac{\partial u}{\partial x}\bigg|_{x_i} \approx \frac{u_i^n - u_{i-1}^n}{\Delta x}$$

The time update is

$$u_i^{n+1}=u_i^n-\frac{a \Delta t}{\Delta x}\left(u_i^n - u_{i-1}^n\right)$$

#### Case $a < 0$

$$\frac{\partial u}{\partial x}\bigg|_{x_i}\approx\frac{u_{i+1}^n - u_i^n}{\Delta x}$$

### Conservative finite volume form

For a conservation law

$$\frac{\partial u}{\partial t} + \frac{\partial f(u)}{\partial x} = 0$$

the upwind numerical flux at interface $i+1/2$ is

$$F_{i+1/2} =
\begin{cases}
f(u_i), & f'(u) > 0 \\
f(u_{i+1}), & f'(u) < 0
\end{cases}$$

For linear advection $f(u) = a u$, this reduces to the previous formulas.

### Stability condition

The upwind scheme is stable under the CFL condition

$$\frac{|a|\Delta t}{\Delta x} \le 1$$

With that, the method is first order accurate in space and in time. However, it introduces numerical diffusion proportional to $\Delta x$. The main advantage is that its unconditionally stable for hyperbolic problems and rather simple to implement. The disadvantage is low accuracy and smoothed gradients at shocks.

For our specific problem of the 1D Vlasov-Poisson equation, we have to split the advection in two steps, one to advance in velocity space and one to advance in position space, where we use different upwind conditions. We can achieve that in python by looping over all position and velocity cells and deciding case per case which sign to choose.


```python
def _advance_velocity(f, dx, dt, nx, nv, v):
    for i in range(nx - 1):
        for j in range(nv - 1):
            if v[j] > 0:
                f[i, j] = f[i, j] - dt / dx * (f[i, j] * v[j] - f[i - 1, j] * v[j])
            else:
                f[i, j] = f[i, j] - dt / dx * (f[i + 1, j] * v[j] - f[i, j] * v[j])
    return f


def _advance_position(f, dv, dt, nx, nv, E):
    for i in range(nx - 1):
        for j in range(nv - 1):
            if E[i] > 0:
                f[i, j] = f[i, j] - dt / dv * (f[i, j] - f[i, j - 1]) * E[i]
            else:
                f[i, j] = f[i, j] - dt / dv * (f[i, j + 1] - f[i, j]) * E[i]
    return f
```

However, this is very slow in python. We can significantly increase the speed by using numpy vectorization. This obscures the readability of the code, however results in a performance boost of approximately 100 times as fast as simple loops. The updates than look like the following functions



```python
def advance_velocity(f, v, dx, dt, v_mask, v_mask_not):
    f[:, v_mask] = f[:, v_mask] - dt / dx * (f - np.roll(f, 1, axis=0))[:, v_mask] * v[v_mask]
    f[:, v_mask_not] = f[:, v_mask_not] - dt / dx * (np.roll(f, -1, axis=0) - f)[:, v_mask_not] * v[v_mask_not]
    return f


def advance_position(f, E, dv, dt, E_mask, E_mask_not):
    f[E_mask, :] = f[E_mask, :] - dt / dv * (f - np.roll(f, 1, axis=1))[E_mask, :] * E[E_mask][..., np.newaxis]
    f[E_mask_not, :] = (
        f[E_mask_not, :] - dt / dv * ((np.roll(f, -1, axis=1) - f)[E_mask_not, :]) * E[E_mask_not][..., np.newaxis]
    )
    return f

```

Finally we can put everything together and implement the time loop. We store the results of every m'th time step into an array for later processing, such as creating a gif of the evolution.



```python
def get_initial_distribution(params):
    # Allocate memory for the distribution function
    f = np.zeros((params["nx"], params["nv"]))
    x = np.linspace(0, params["L"], params["nx"], endpoint=False)
    v = np.linspace(params["vmin"], params["vmax"], params["nv"])
    dv = np.abs(v[1] - v[0])  # velocity step size
    dx = np.abs(x[1] - x[0])  # spatial step size
    dt = params["CFL"] * min(dx / params["vmax"], dv / 1.0)  # assuming max|E| ~ 1

    print("Simulation parameters:")
    print(f"  Spatial step size (dx): {dx}")
    print(f"  Velocity step size (dv): {dv}")
    print(f"  Time step size (dt): {dt}")

    return f, x, v, dx, dv, dt


def plot_distribution(f, x, v):
    # plt.figure(figsize=(6, 4), dpi=100)
    plt.imshow(f.T, cmap=CMAP, aspect="auto", origin="lower", extent=[x.min(), x.max(), v.min(), v.max()])
    # plt.colorbar()
    plt.title("Distribution function $f(x,v)$")
    plt.xlabel("Position $x$")
    plt.ylabel("Velocity $v$")
    # plt.show()
    return


def run_simulation(f, x, v, dx, dv, dt, params):
    # Allocate memory for the results
    result = np.zeros((params["maxiter"] // params["save_interval"], params["nx"], params["nv"]))
    E_k_hist = []

    # Get masks for upwind schemes
    v_mask = np.argwhere(v > 0)
    v_mask_not = np.argwhere(v <= 0)

    for it in tqdm.tqdm(range(params["maxiter"])):
        # Advance velocity distribution
        f_new = advance_velocity(f, v, dx, dt, v_mask, v_mask_not)

        # Integrate f along velocity axis to get charge density
        rho = params["q_e"] * np.trapezoid(f_new, dx=dv, axis=1)

        # Add constant ion background
        rho = rho - np.mean(rho)

        # Solve Poisson equation
        E = solve_poisson_fft(rho, dx, params["nx"])

        # Get mask for upwind scheme
        E_mask = np.argwhere(E > 0)
        E_mask_not = np.argwhere(E <= 0)

        # Advance position distribution
        f_new = advance_position(f_new, E, dv, dt, E_mask, E_mask_not)

        if params["export_mode"]:
            # Solve Poisson with FFT
            rho_k = np.fft.fft(rho)
            k_vals = 2 * np.pi * np.fft.fftfreq(params["nx"], d=dx)

            # Avoid division by zero at k=0
            E_k = np.zeros_like(rho_k, dtype=complex)
            nonzero = k_vals != 0
            E_k[nonzero] = (1j / k_vals[nonzero]) * rho_k[nonzero]
            E = np.real(np.fft.ifft(E_k))

            # Save fundamental mode (k=±1 depending on domain length L)
            # If domain length is L=2π/k0, then fundamental mode index is 1
            E1 = E_k[1] / params["nx"]
            E_k_hist.append(E1)

        # Save every mth step
        if it % params["save_interval"] == 0:
            result[it // params["save_interval"]] = f.copy()
            if params["plot"]:
                plot_distribution(f, x, v)
                plt.savefig(f"/home/max/Temp/{params['name']}_{it:05d}.png")
    return result, E_k_hist
```

## Kinetic instabilities in collisionless plasmas

In collisionless plasmas, wave propagation and stability are described by the Vlasov equation coupled to Poisson’s equation. Landau damping, the two stream instability, and the bump on tail instability arise from resonant wave particle interactions.

Linearizing around a homogeneous equilibrium $f_0(v)$ leads to the electrostatic dispersion relation

$$1 + \frac{1}{k^2 \varepsilon_0}\sum_s \frac{q_s^2}{m_s}\int \frac{\partial f_{0s}/\partial v}{v - \omega/k} \, dv= 0$$

## Landau damping

Landau damping occurs when the equilibrium distribution function $f_0(v)$ is monotonic decreasing. The imaginary part of the wave frequency is obtained by evaluating the dispersion relation using the Landau contour, yielding

$$\gamma =\operatorname{Im}(\omega)=-\frac{\pi}{2}\frac{\omega_p^2}{k}\left.\frac{\partial f_0}{\partial v}\right|_{v = \omega_r/k}$$

where $\omega_p$ is the plasma frequency and $\omega_r$ is the real part of the wave frequency. For $\partial f_0 / \partial v < 0$, the damping rate $\gamma$ is negative, corresponding to exponential decay of the wave amplitude.



```python
params = {
    "nx": 128 + 1,  # Number of spatial positions
    "nv": 128 + 1,  # Number of velocity values
    "vmin": -10,  # Minimal possible velocity
    "vmax": 10,  # Maximal possible velocity
    "q_e": -1,  # Electron charge
    "n_particles": 1,  # Initial particle number
    "k": 0.1,  # wave number of perturbation
    "L": 2 * np.pi / 0.1,  # domain length for k
    "CFL": 0.5,  # CFL number
    "maxiter": 4000,  # Maximum number of iterations
    "plot": True,  # Whether to plot intermediate results
    "save_interval": 100,  # Save every mth step
    "export_mode": True,  # Whether to export mode data
    "save_path": "/home/max/Temp",  # Path to save output files
    "name": "landau_damping",  # Simulation name
}


def maxwellian(v, vbeam, vth):
    return np.exp(-0.5 * ((v - vbeam) / vth) ** 2) / (np.sqrt(2 * np.pi) * vth)


# Get initial distribution and simulation parameters
f, x, v, dx, dv, dt = get_initial_distribution(params)

# Maxwellian
vth = 3.0
f0 = maxwellian(v, 0.0, vth)

# Small perturbation
alpha = 1e-1
f = np.zeros((params["nx"], params["nv"]))
for i in range(params["nx"]):
    f[i, :] = f0 * (1 + alpha * np.cos(params["k"] * x[i]))

# Normalize the distribution
f = f / f.sum()

plot_distribution(f, x, v)
```

    Simulation parameters:
      Spatial step size (dx): 0.4870686284635338
      Velocity step size (dv): 0.15625
      Time step size (dt): 0.02435343142317669



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_12_1.png)
    



```python
result, E_k_hist = run_simulation(f, x, v, dx, dv, dt, params)
```

      0%|          | 0/4000 [00:00<?, ?it/s]

    100%|██████████| 4000/4000 [00:38<00:00, 102.75it/s]



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_13_2.png)
    



```python
E_k_hist = np.array(E_k_hist)
t = np.arange(len(E_k_hist)) * dt

# Amplitude of first Fourier mode
amp = np.abs(E_k_hist)

plt.plot(t, amp)
plt.yscale("log")
plt.xlabel("Time")
plt.ylabel(r"$|E_{k=0.5}(t)|$")
plt.title("Landau damping: electric field amplitude")
plt.show()

```


    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_14_0.png)
    



```python
# find peaks
peaks, _ = scipy.signal.find_peaks(amp)
plt.plot(t[peaks], amp[peaks], "x")
plt.yscale("log")
plt.xscale("log")
plt.xlabel("Time")
plt.ylabel(r"$|E_k(t)|$")
plt.grid(which="both", linestyle="--", linewidth=0.5)
plt.show()

# get slope
fit_range = (t > 0.5) & (t < 30)
slope, intercept, *_ = scipy.stats.linregress(t[peaks], np.log(amp[peaks]))

print(f"Numerical damping rate gamma = {slope:.4f}")
```


    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_15_0.png)
    


    Numerical damping rate gamma = -0.0184


## Two stream instability

The two stream instability arises when the equilibrium distribution consists of two drifting populations. A simple cold plasma model uses

$$f_0(v) =\frac{n_0}{2}\left[\delta(v - v_0) + \delta(v + v_0)\right]$$

The resulting dispersion relation is

$$1 - \frac{\omega_p^2}{(\omega - k v_0)^2}- \frac{\omega_p^2}{(\omega + k v_0)^2}= 0$$

For sufficiently large drift velocity $v_0$, this equation admits solutions with $\operatorname{Im}(\omega) > 0$, indicating exponential growth of electric field perturbations.



```python
params = {
    "nx": 256 + 1,
    "nv": 256 + 1,
    "vmin": -10,
    "vmax": 10,
    "q_e": -1000,
    "n_particles": 1,
    "k": 0.1,
    "L": 256,
    "CFL": 0.5,
    "maxiter": 4000,
    "plot": True,
    "save_interval": 100,
    "export_mode": True,
    "save_path": "/home/max/Temp",
    "name": "two_stream_instability",
}

# Get initial distribution and simulation parameters
f, x, v, dx, dv, dt = get_initial_distribution(params)

vth = 0.01 * params["vmax"]  # Thermal velocity
vbeam = 0.2 * params["vmax"]  # Beam velocity

# Initial distribution for two-stream instability
for i in range(params["nv"]):
    f[:, i] = np.exp(-((v[i] - vbeam) ** 2) / vth**2) + np.exp(-((v[i] + vbeam) ** 2) / vth**2)

# Add noise
f += np.random.normal(0, 0.01, f.shape)

# Normalize the distribution
f = f / f.sum()

# Plot initial distribution
plot_distribution(f, x, v)
```

    Simulation parameters:
      Spatial step size (dx): 0.9961089494163424
      Velocity step size (dv): 0.078125
      Time step size (dt): 0.0390625



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_17_1.png)
    



```python
result, E_k_hist = run_simulation(f, x, v, dx, dv, dt, params)
```

      0%|          | 0/4000 [00:00<?, ?it/s]

    100%|██████████| 4000/4000 [00:43<00:00, 91.53it/s] 



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_18_2.png)
    


## Bump on tail instability

The bump on tail instability occurs when the distribution function has a local maximum at high velocity, so that

$$\left.\frac{\partial f_0}{\partial v}\right|_{v = \omega/k}> 0$$

The growth rate follows the same expression as Landau damping but with opposite sign,

$$\gamma =-\frac{\pi}{2}\frac{\omega_p^2}{k}\left.\frac{\partial f_0}{\partial v}\right|_{v = \omega_r/k}$$

A positive slope at the resonant velocity leads to $\gamma > 0$, corresponding to wave amplification.



```python
params = {
    "nx": 256 + 1,
    "nv": 256 + 1,
    "vmin": -10,
    "vmax": 10,
    "q_e": -1,
    "n_particles": 1,
    "k": 0.1,
    "L": 2 * np.pi / 0.1,
    "CFL": 0.5,
    "maxiter": 4000,
    "plot": True,
    "save_interval": 100,
    "export_mode": True,
    "save_path": "/home/max/Temp",
    "name": "bump_on_tail",
}

# Get initial distribution and simulation parameters
f, x, v, dx, dv, dt = get_initial_distribution(params)

vth = 0.2 * params["vmax"]  # Thermal velocity
vbeam = 0.6 * params["vmax"]  # Beam velocity

alpha_b = 0.1
eps = 1e-4

fM = lambda v, u, vt: np.exp(-((v - u) ** 2) / (2 * vt**2)) / (np.sqrt(2 * np.pi) * vt)

f0v = (1 - alpha_b) * fM(v, 0, vth) + alpha_b * fM(v, vbeam, 0.3)

f = np.zeros((params["nx"], params["nv"]))
for i in range(params["nx"]):
    f[i, :] = f0v * (1 + eps * np.cos(10 * x[i]))

f += np.random.normal(0, 0.001, f.shape)

f /= np.sum(f[0, :]) * dv

# Plot initial distribution
plot_distribution(f, x, v)
```

    Simulation parameters:
      Spatial step size (dx): 0.24448191856729906
      Velocity step size (dv): 0.078125
      Time step size (dt): 0.012224095928364953



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_20_1.png)
    



```python
result, E_k_hist = run_simulation(f, x, v, dx, dv, dt, params)
```

      0%|          | 0/4000 [00:00<?, ?it/s]

    100%|██████████| 4000/4000 [00:43<00:00, 92.12it/s] 



    
![png](vlasov_simulation_fft_files/vlasov_simulation_fft_21_2.png)
    


### Summary

All three phenomena are described by the same linearized Vlasov dispersion relation. Landau damping corresponds to negative slope of the distribution function at resonance, while the two stream and bump on tail instabilities arise from nonmonotonic velocity distributions that provide free energy for wave growth.

The jupyter notebook of this tutorial can be downloaded ![here](/docs/vlasov_simulation_fft.ipynb)



```python
import subprocess

def create_video(name):
    subprocess.run(["ffmpeg", "-framerate", "15", "-pattern_type", "glob", "-i", f"/home/max/Temp/{name}_*.png", "-c:v", "libx264", "-preset", "slow", "-crf", "18", "-pix_fmt", "yuv420p", f"/home/max/Repos/portfolio/static/videos/{name}.mp4"])

create_video("landau_damping")
create_video("two_stream_instability")
create_video("bump_on_tail")
```

    ffmpeg version n8.0.1 Copyright (c) 2000-2025 the FFmpeg developers
      built with gcc 15.2.1 (GCC) 20251112
      configuration: --prefix=/usr --disable-debug --disable-static --disable-stripping --enable-amf --enable-avisynth --enable-cuda-llvm --enable-lto --enable-fontconfig --enable-frei0r --enable-gmp --enable-gnutls --enable-gpl --enable-ladspa --enable-libaom --enable-libass --enable-libbluray --enable-libbs2b --enable-libdav1d --enable-libdrm --enable-libdvdnav --enable-libdvdread --enable-libfreetype --enable-libfribidi --enable-libglslang --enable-libgsm --enable-libharfbuzz --enable-libiec61883 --enable-libjack --enable-libjxl --enable-libmodplug --enable-libmp3lame --enable-libopencore_amrnb --enable-libopencore_amrwb --enable-libopenjpeg --enable-libopenmpt --enable-libopus --enable-libplacebo --enable-libpulse --enable-librav1e --enable-librsvg --enable-librubberband --enable-libsnappy --enable-libsoxr --enable-libspeex --enable-libsrt --enable-libssh --enable-libsvtav1 --enable-libtheora --enable-libv4l2 --enable-libvidstab --enable-libvmaf --enable-libvorbis --enable-libvpl --enable-libvpx --enable-libwebp --enable-libx264 --enable-libx265 --enable-libxcb --enable-libxml2 --enable-libxvid --enable-libzimg --enable-libzmq --enable-nvdec --enable-nvenc --enable-opencl --enable-opengl --enable-shared --enable-vapoursynth --enable-version3 --enable-vulkan
      libavutil      60.  8.100 / 60.  8.100
      libavcodec     62. 11.100 / 62. 11.100
      libavformat    62.  3.100 / 62.  3.100
      libavdevice    62.  1.100 / 62.  1.100
      libavfilter    11.  4.100 / 11.  4.100
      libswscale      9.  1.100 /  9.  1.100
      libswresample   6.  1.100 /  6.  1.100
    Input #0, image2, from '/home/max/Temp/landau_damping_*.png':
      Duration: 00:00:02.67, start: 0.000000, bitrate: N/A
      Stream #0:0: Video: png, rgba(pc, gbr/unknown/unknown), 600x400 [SAR 3937:3937 DAR 3:2], 15 fps, 15 tbr, 15 tbn
    Stream mapping:
      Stream #0:0 -> #0:0 (png (native) -> h264 (libx264))
    Press [q] to stop, [?] for help
    [libx264 @ 0x5629c40d49c0] using SAR=1/1
    [libx264 @ 0x5629c40d49c0] using cpu capabilities: MMX2 SSE2Fast SSSE3 SSE4.2 AVX
    [libx264 @ 0x5629c40d49c0] profile High, level 2.2, 4:2:0, 8-bit
    [libx264 @ 0x5629c40d49c0] 264 - core 165 r3222 b35605a - H.264/MPEG-4 AVC codec - Copyleft 2003-2025 - http://www.videolan.org/x264.html - options: cabac=1 ref=5 deblock=1:0:0 analyse=0x3:0x113 me=hex subme=8 psy=1 psy_rd=1.00:0.00 mixed_ref=1 me_range=16 chroma_me=1 trellis=2 8x8dct=1 cqm=0 deadzone=21,11 fast_pskip=1 chroma_qp_offset=-2 threads=12 lookahead_threads=2 sliced_threads=0 nr=0 decimate=1 interlaced=0 bluray_compat=0 constrained_intra=0 bframes=3 b_pyramid=2 b_adapt=1 b_bias=0 direct=3 weightb=1 open_gop=0 weightp=2 keyint=250 keyint_min=15 scenecut=40 intra_refresh=0 rc_lookahead=50 rc=crf mbtree=1 crf=18.0 qcomp=0.60 qpmin=0 qpmax=69 qpstep=4 ip_ratio=1.40 aq=1:1.00
    Output #0, mp4, to '/home/max/Repos/portfolio/static/videos/landau_damping.mp4':
      Metadata:
        encoder         : Lavf62.3.100
      Stream #0:0: Video: h264 (avc1 / 0x31637661), yuv420p(tv, progressive), 600x400 [SAR 1:1 DAR 3:2], q=2-31, 15 fps, 15360 tbn
        Metadata:
          encoder         : Lavc62.11.100 libx264
        Side data:
          cpb: bitrate max/min/avg: 0/0/0 buffer size: 0 vbv_delay: N/A
    [out#0/mp4 @ 0x5629c40d38c0] video:85KiB audio:0KiB subtitle:0KiB other streams:0KiB global headers:0KiB muxing overhead: 1.551387%
    frame=   40 fps=0.0 q=-1.0 Lsize=      86KiB time=00:00:02.53 bitrate= 278.6kbits/s speed= 4.4x elapsed=0:00:00.57    
    [libx264 @ 0x5629c40d49c0] frame I:1     Avg QP:16.57  size: 14402
    [libx264 @ 0x5629c40d49c0] frame P:10    Avg QP:17.07  size:  3343
    [libx264 @ 0x5629c40d49c0] frame B:29    Avg QP:18.35  size:  1323
    [libx264 @ 0x5629c40d49c0] consecutive B-frames:  2.5%  0.0%  7.5% 90.0%
    [libx264 @ 0x5629c40d49c0] mb I  I16..4: 52.2% 32.9% 14.8%
    [libx264 @ 0x5629c40d49c0] mb P  I16..4: 10.0%  6.0%  3.0%  P16..4: 12.0% 12.8%  2.1%  0.0%  0.0%    skip:54.1%
    [libx264 @ 0x5629c40d49c0] mb B  I16..4:  1.7%  0.6%  0.0%  B16..8: 27.4%  5.5%  0.4%  direct: 3.4%  skip:60.9%  L0:46.7% L1:31.4% BI:22.0%
    [libx264 @ 0x5629c40d49c0] 8x8 transform intra:31.2% inter:82.2%
    [libx264 @ 0x5629c40d49c0] direct mvs  spatial:79.3% temporal:20.7%
    [libx264 @ 0x5629c40d49c0] coded y,uvDC,uvAC intra: 33.8% 73.0% 53.1% inter: 5.2% 21.3% 5.0%
    [libx264 @ 0x5629c40d49c0] i16 v,h,dc,p:  5% 85%  2%  8%
    [libx264 @ 0x5629c40d49c0] i8 v,h,dc,ddl,ddr,vr,hd,vl,hu: 11% 74% 11%  0%  0%  0%  1%  0%  2%
    [libx264 @ 0x5629c40d49c0] i4 v,h,dc,ddl,ddr,vr,hd,vl,hu: 13% 71%  9%  1%  1%  1%  1%  1%  3%
    [libx264 @ 0x5629c40d49c0] i8c dc,h,v,p:  9% 88%  2%  2%
    [libx264 @ 0x5629c40d49c0] Weighted P-Frames: Y:0.0% UV:0.0%
    [libx264 @ 0x5629c40d49c0] ref P L0: 56.1%  6.2% 23.0%  5.1%  9.7%
    [libx264 @ 0x5629c40d49c0] ref B L0: 69.0% 19.0% 10.8%  1.2%
    [libx264 @ 0x5629c40d49c0] ref B L1: 95.0%  5.0%
    [libx264 @ 0x5629c40d49c0] kb/s:258.60
    ffmpeg version n8.0.1 Copyright (c) 2000-2025 the FFmpeg developers
      built with gcc 15.2.1 (GCC) 20251112
      configuration: --prefix=/usr --disable-debug --disable-static --disable-stripping --enable-amf --enable-avisynth --enable-cuda-llvm --enable-lto --enable-fontconfig --enable-frei0r --enable-gmp --enable-gnutls --enable-gpl --enable-ladspa --enable-libaom --enable-libass --enable-libbluray --enable-libbs2b --enable-libdav1d --enable-libdrm --enable-libdvdnav --enable-libdvdread --enable-libfreetype --enable-libfribidi --enable-libglslang --enable-libgsm --enable-libharfbuzz --enable-libiec61883 --enable-libjack --enable-libjxl --enable-libmodplug --enable-libmp3lame --enable-libopencore_amrnb --enable-libopencore_amrwb --enable-libopenjpeg --enable-libopenmpt --enable-libopus --enable-libplacebo --enable-libpulse --enable-librav1e --enable-librsvg --enable-librubberband --enable-libsnappy --enable-libsoxr --enable-libspeex --enable-libsrt --enable-libssh --enable-libsvtav1 --enable-libtheora --enable-libv4l2 --enable-libvidstab --enable-libvmaf --enable-libvorbis --enable-libvpl --enable-libvpx --enable-libwebp --enable-libx264 --enable-libx265 --enable-libxcb --enable-libxml2 --enable-libxvid --enable-libzimg --enable-libzmq --enable-nvdec --enable-nvenc --enable-opencl --enable-opengl --enable-shared --enable-vapoursynth --enable-version3 --enable-vulkan
      libavutil      60.  8.100 / 60.  8.100
      libavcodec     62. 11.100 / 62. 11.100
      libavformat    62.  3.100 / 62.  3.100
      libavdevice    62.  1.100 / 62.  1.100
      libavfilter    11.  4.100 / 11.  4.100
      libswscale      9.  1.100 /  9.  1.100
      libswresample   6.  1.100 /  6.  1.100
    Input #0, image2, from '/home/max/Temp/two_stream_instability_*.png':
      Duration: 00:00:02.67, start: 0.000000, bitrate: N/A
      Stream #0:0: Video: png, rgba(pc, gbr/unknown/unknown), 600x400 [SAR 3937:3937 DAR 3:2], 15 fps, 15 tbr, 15 tbn
    Stream mapping:
      Stream #0:0 -> #0:0 (png (native) -> h264 (libx264))
    Press [q] to stop, [?] for help
    [libx264 @ 0x55d4dafb7e40] using SAR=1/1
    [libx264 @ 0x55d4dafb7e40] using cpu capabilities: MMX2 SSE2Fast SSSE3 SSE4.2 AVX
    [libx264 @ 0x55d4dafb7e40] profile High, level 2.2, 4:2:0, 8-bit
    [libx264 @ 0x55d4dafb7e40] 264 - core 165 r3222 b35605a - H.264/MPEG-4 AVC codec - Copyleft 2003-2025 - http://www.videolan.org/x264.html - options: cabac=1 ref=5 deblock=1:0:0 analyse=0x3:0x113 me=hex subme=8 psy=1 psy_rd=1.00:0.00 mixed_ref=1 me_range=16 chroma_me=1 trellis=2 8x8dct=1 cqm=0 deadzone=21,11 fast_pskip=1 chroma_qp_offset=-2 threads=12 lookahead_threads=2 sliced_threads=0 nr=0 decimate=1 interlaced=0 bluray_compat=0 constrained_intra=0 bframes=3 b_pyramid=2 b_adapt=1 b_bias=0 direct=3 weightb=1 open_gop=0 weightp=2 keyint=250 keyint_min=15 scenecut=40 intra_refresh=0 rc_lookahead=50 rc=crf mbtree=1 crf=18.0 qcomp=0.60 qpmin=0 qpmax=69 qpstep=4 ip_ratio=1.40 aq=1:1.00
    Output #0, mp4, to '/home/max/Repos/portfolio/static/videos/two_stream_instability.mp4':
      Metadata:
        encoder         : Lavf62.3.100
      Stream #0:0: Video: h264 (avc1 / 0x31637661), yuv420p(tv, progressive), 600x400 [SAR 1:1 DAR 3:2], q=2-31, 15 fps, 15360 tbn
        Metadata:
          encoder         : Lavc62.11.100 libx264
        Side data:
          cpb: bitrate max/min/avg: 0/0/0 buffer size: 0 vbv_delay: N/A
    [out#0/mp4 @ 0x55d4dafb72c0] video:181KiB audio:0KiB subtitle:0KiB other streams:0KiB global headers:0KiB muxing overhead: 0.730619%
    frame=   40 fps=0.0 q=-1.0 Lsize=     183KiB time=00:00:02.53 bitrate= 590.4kbits/s speed=4.69x elapsed=0:00:00.53    
    [libx264 @ 0x55d4dafb7e40] frame I:1     Avg QP:16.02  size: 15619
    [libx264 @ 0x55d4dafb7e40] frame P:11    Avg QP:15.93  size:  6634
    [libx264 @ 0x55d4dafb7e40] frame B:28    Avg QP:16.55  size:  3440
    [libx264 @ 0x55d4dafb7e40] consecutive B-frames:  5.0%  5.0%  0.0% 90.0%
    [libx264 @ 0x55d4dafb7e40] mb I  I16..4: 41.3% 45.1% 13.7%
    [libx264 @ 0x55d4dafb7e40] mb P  I16..4:  4.9%  9.2%  2.2%  P16..4:  8.0%  7.8%  4.6%  0.0%  0.0%    skip:63.3%
    [libx264 @ 0x55d4dafb7e40] mb B  I16..4:  0.4%  1.5%  0.3%  B16..8: 15.0%  7.9%  2.5%  direct: 3.8%  skip:68.5%  L0:34.7% L1:33.4% BI:32.0%
    [libx264 @ 0x55d4dafb7e40] 8x8 transform intra:55.1% inter:88.6%
    [libx264 @ 0x55d4dafb7e40] direct mvs  spatial:75.0% temporal:25.0%
    [libx264 @ 0x55d4dafb7e40] coded y,uvDC,uvAC intra: 58.3% 69.8% 60.7% inter: 11.5% 22.1% 13.3%
    [libx264 @ 0x55d4dafb7e40] i16 v,h,dc,p: 34% 57%  5%  4%
    [libx264 @ 0x55d4dafb7e40] i8 v,h,dc,ddl,ddr,vr,hd,vl,hu:  9% 30% 17%  7%  6%  4% 10%  5% 12%
    [libx264 @ 0x55d4dafb7e40] i4 v,h,dc,ddl,ddr,vr,hd,vl,hu: 19% 36% 17%  3%  5%  3%  6%  4%  7%
    [libx264 @ 0x55d4dafb7e40] i8c dc,h,v,p: 26% 56% 13%  5%
    [libx264 @ 0x55d4dafb7e40] Weighted P-Frames: Y:0.0% UV:0.0%
    [libx264 @ 0x55d4dafb7e40] ref P L0: 56.3%  7.4% 19.8%  9.8%  6.7%
    [libx264 @ 0x55d4dafb7e40] ref B L0: 85.7% 10.6%  3.0%  0.7%
    [libx264 @ 0x55d4dafb7e40] ref B L1: 97.0%  3.0%
    [libx264 @ 0x55d4dafb7e40] kb/s:554.72
    ffmpeg version n8.0.1 Copyright (c) 2000-2025 the FFmpeg developers
      built with gcc 15.2.1 (GCC) 20251112
      configuration: --prefix=/usr --disable-debug --disable-static --disable-stripping --enable-amf --enable-avisynth --enable-cuda-llvm --enable-lto --enable-fontconfig --enable-frei0r --enable-gmp --enable-gnutls --enable-gpl --enable-ladspa --enable-libaom --enable-libass --enable-libbluray --enable-libbs2b --enable-libdav1d --enable-libdrm --enable-libdvdnav --enable-libdvdread --enable-libfreetype --enable-libfribidi --enable-libglslang --enable-libgsm --enable-libharfbuzz --enable-libiec61883 --enable-libjack --enable-libjxl --enable-libmodplug --enable-libmp3lame --enable-libopencore_amrnb --enable-libopencore_amrwb --enable-libopenjpeg --enable-libopenmpt --enable-libopus --enable-libplacebo --enable-libpulse --enable-librav1e --enable-librsvg --enable-librubberband --enable-libsnappy --enable-libsoxr --enable-libspeex --enable-libsrt --enable-libssh --enable-libsvtav1 --enable-libtheora --enable-libv4l2 --enable-libvidstab --enable-libvmaf --enable-libvorbis --enable-libvpl --enable-libvpx --enable-libwebp --enable-libx264 --enable-libx265 --enable-libxcb --enable-libxml2 --enable-libxvid --enable-libzimg --enable-libzmq --enable-nvdec --enable-nvenc --enable-opencl --enable-opengl --enable-shared --enable-vapoursynth --enable-version3 --enable-vulkan
      libavutil      60.  8.100 / 60.  8.100
      libavcodec     62. 11.100 / 62. 11.100
      libavformat    62.  3.100 / 62.  3.100
      libavdevice    62.  1.100 / 62.  1.100
      libavfilter    11.  4.100 / 11.  4.100
      libswscale      9.  1.100 /  9.  1.100
      libswresample   6.  1.100 /  6.  1.100
    Input #0, image2, from '/home/max/Temp/bump_on_tail_*.png':
      Duration: 00:00:02.67, start: 0.000000, bitrate: N/A
      Stream #0:0: Video: png, rgba(pc, gbr/unknown/unknown), 600x400 [SAR 3937:3937 DAR 3:2], 15 fps, 15 tbr, 15 tbn
    Stream mapping:
      Stream #0:0 -> #0:0 (png (native) -> h264 (libx264))
    Press [q] to stop, [?] for help
    [libx264 @ 0x55d93e024d00] using SAR=1/1
    [libx264 @ 0x55d93e024d00] using cpu capabilities: MMX2 SSE2Fast SSSE3 SSE4.2 AVX
    [libx264 @ 0x55d93e024d00] profile High, level 2.2, 4:2:0, 8-bit
    [libx264 @ 0x55d93e024d00] 264 - core 165 r3222 b35605a - H.264/MPEG-4 AVC codec - Copyleft 2003-2025 - http://www.videolan.org/x264.html - options: cabac=1 ref=5 deblock=1:0:0 analyse=0x3:0x113 me=hex subme=8 psy=1 psy_rd=1.00:0.00 mixed_ref=1 me_range=16 chroma_me=1 trellis=2 8x8dct=1 cqm=0 deadzone=21,11 fast_pskip=1 chroma_qp_offset=-2 threads=12 lookahead_threads=2 sliced_threads=0 nr=0 decimate=1 interlaced=0 bluray_compat=0 constrained_intra=0 bframes=3 b_pyramid=2 b_adapt=1 b_bias=0 direct=3 weightb=1 open_gop=0 weightp=2 keyint=250 keyint_min=15 scenecut=40 intra_refresh=0 rc_lookahead=50 rc=crf mbtree=1 crf=18.0 qcomp=0.60 qpmin=0 qpmax=69 qpstep=4 ip_ratio=1.40 aq=1:1.00
    Output #0, mp4, to '/home/max/Repos/portfolio/static/videos/bump_on_tail.mp4':
      Metadata:
        encoder         : Lavf62.3.100
      Stream #0:0: Video: h264 (avc1 / 0x31637661), yuv420p(tv, progressive), 600x400 [SAR 1:1 DAR 3:2], q=2-31, 15 fps, 15360 tbn
        Metadata:
          encoder         : Lavc62.11.100 libx264
        Side data:
          cpb: bitrate max/min/avg: 0/0/0 buffer size: 0 vbv_delay: N/A
    [out#0/mp4 @ 0x55d93e023cc0] video:56KiB audio:0KiB subtitle:0KiB other streams:0KiB global headers:0KiB muxing overhead: 2.342189%
    frame=   40 fps=0.0 q=-1.0 Lsize=      58KiB time=00:00:02.53 bitrate= 186.0kbits/s speed=8.11x elapsed=0:00:00.31    
    [libx264 @ 0x55d93e024d00] frame I:1     Avg QP:16.67  size: 12953
    [libx264 @ 0x55d93e024d00] frame P:12    Avg QP:16.41  size:  2159
    [libx264 @ 0x55d93e024d00] frame B:27    Avg QP:18.71  size:   667
    [libx264 @ 0x55d93e024d00] consecutive B-frames: 10.0%  0.0%  0.0% 90.0%
    [libx264 @ 0x55d93e024d00] mb I  I16..4: 55.4% 27.7% 16.9%
    [libx264 @ 0x55d93e024d00] mb P  I16..4:  3.5%  1.8%  0.2%  P16..4: 13.8%  8.8%  4.6%  0.0%  0.0%    skip:67.3%
    [libx264 @ 0x55d93e024d00] mb B  I16..4:  0.7%  0.1%  0.0%  B16..8: 17.3%  3.0%  0.4%  direct: 1.2%  skip:77.3%  L0:43.4% L1:44.1% BI:12.5%
    [libx264 @ 0x55d93e024d00] 8x8 transform intra:27.7% inter:67.4%
    [libx264 @ 0x55d93e024d00] direct mvs  spatial:77.8% temporal:22.2%
    [libx264 @ 0x55d93e024d00] coded y,uvDC,uvAC intra: 28.2% 39.5% 26.8% inter: 3.4% 10.6% 2.7%
    [libx264 @ 0x55d93e024d00] i16 v,h,dc,p: 11% 82%  4%  3%
    [libx264 @ 0x55d93e024d00] i8 v,h,dc,ddl,ddr,vr,hd,vl,hu: 19% 40% 23%  1%  2%  1%  5%  1%  7%
    [libx264 @ 0x55d93e024d00] i4 v,h,dc,ddl,ddr,vr,hd,vl,hu: 27% 31% 21%  2%  3%  3%  4%  3%  7%
    [libx264 @ 0x55d93e024d00] i8c dc,h,v,p: 16% 78%  4%  1%
    [libx264 @ 0x55d93e024d00] Weighted P-Frames: Y:0.0% UV:0.0%
    [libx264 @ 0x55d93e024d00] ref P L0: 61.2%  7.1% 21.0%  5.0%  4.7%  1.0%
    [libx264 @ 0x55d93e024d00] ref B L0: 70.2% 16.5% 11.9%  1.3%
    [libx264 @ 0x55d93e024d00] ref B L1: 95.4%  4.6%
    [libx264 @ 0x55d93e024d00] kb/s:170.59



```python

```
