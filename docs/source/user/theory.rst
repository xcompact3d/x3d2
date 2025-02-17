Theoretical Background
======================

.. _theoretical-background:

Numerical method framework
--------------------------

.. _numerical_method:

Governing equations
~~~~~~~~~~~~~~~~~~~

.. _governing-equations:

The governing equations are the forced incompressible Navier-Stokes equations:


.. math::
   :label: incomp-ns

   \frac{\partial(\mathbf{u})}{\partial{t}} = -\nabla{p} - \frac{1}{2}[\nabla (\mathbf{u} \otimes \mathbf{u})
   + (\mathbf{u} \cdot \nabla) \mathbf{u}] + \nu \nabla^2\mathbf{u} + \mathbf{f}

.. math::
   :label: div-free

    \nabla \cdot \mathbf{u} = 0

where :math:`\mathbf{u}` is the velocity field, :math:`p` is the pressure field, :math:`\nu` is the kinematic viscosity,
and :math:`\mathbf{f}` is the external force field. Eq. :eq:`incomp-ns` is the momentum equation and  Eq. :eq:`div-free` 
is the incompressibility constraint. Further details on the numerical methods can be found in :cite:`bartholomew_sx_20, laizet_jcp_09, lamballais_jcp_11`.

Time advancement
~~~~~~~~~~~~~~~~

The time advancement of Eq. :eq:`incomp-ns` can be expressed as:

.. math::
    :label: time-adv

    &\frac{\mathbf{u}^* - \mathbf{u}^k}{\Delta{t}} = a_k\mathbf{F}^k + b_k\mathbf{F}^{k-1} - c_k\nabla\tilde{p}^k + c_k\tilde{\mathbf{f}}^{k+1} \\
    &\frac{\mathbf{u}^{**} - \mathbf{u}^*}{\Delta{t}} = c_k\nabla\tilde{p}^k\\
    &\frac{\mathbf{u}^{k+1} - \mathbf{u}^{**}}{\Delta{t}} = -c_k\nabla\tilde{p}^{k+1}

with 

.. math::
   :label: incomp-ns-rhs

   \mathbf{F}^k = - \frac{1}{2}[\nabla (\mathbf{u}^k \otimes \mathbf{u}^k)
   + (\mathbf{u}^k \cdot \nabla) \mathbf{u}^k] + \nu \nabla^2\mathbf{u}^k

and

.. math::
    :label: pressure-correction

    \tilde{p}^{k+1} = \frac{1}{c_k\Delta{t}}\int_{t_k}^{t_{k+1}} p\, \mathrm{d}t, \quad \tilde{\mathbf{f}}^{k+1} = \frac{1}{c_k\Delta{t}}\int_{t_k}^{t_{k+1}} \mathbf{f}\, \mathrm{d}t


for a Runge-Kutta scheme with coefficients :math:`a_k`, :math:`b_k`, and :math:`c_k=a_k+b_k` and :math:`n_k` sub-time steps :math:`k=1,\dots{n_k}` with :math:`t_1=t_n` and :math:`t_{n_k} = t_{n+1}`. Pressure and forcing terms are expressed through their time-averaged values on a given sub-step :math:`c_k\Delta{t}`, indicated by the tilde in :math:`\tilde{p}^{k+1}` and :math:`\tilde{\mathbf{f}}^{k+1}`.

Boundary conditions
~~~~~~~~~~~~~~~~~~~

The governing equations Eq. :eq:`incomp-ns` and :eq:`div-free` are solved in a computational domain :math:`L_x \times L_y \times L_z` discretised on a Cartesian mesh of :math:`n_x \times n_y \times n_z` nodes. 

At the boundaries of the time periodic, free-slip, no-slip, or open conditions can be applied depending on the flow configuration considered. Period and free-slip boundary conditions can be imposed directly via the spatial discretisation without specific care in time advancements. 

However, the use of Dirichlet conditions on the velocity (for no-slip or open conditions) needs to be defined according to time advancement procedure. Conventional homogeneous Neumann conditions are used to solve the pressure.

Role of intermediate velocities
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this approach, we introduce two intermediate velocities :math:`\mathbf{u}^*` and :math:`\mathbf{u}^{**}`. The motivation for using these intermediate steps is to enforce the divergence-free condition at the walls while also satisfying Dirichlet boundary conditions.

First, we compute a velocity field :math:`\mathbf{u}^*` that satisfies the momentum equation without yet enforcing incompressibility. This provides a preliminary estimate of the velocity. Next, we modify :math:`\mathbf{u}^*` by incorporating the pressure gradient from the previous time step :math:`\nabla{p}^k` to obtain :math:`\mathbf{u}^{**}`.

.. math::
    :label: first_vel

    \mathbf{u}^{**}\big|_w = \mathbf{u}^*\big|_w + \Delta{t}\cdot{c_k} \nabla{p}^k

To ensure :math:`\nabla\cdot\mathbf{u}^{k+1}=0` we use the pressure gradient :math:`\nabla{p}^{k+1}` from the current time-step:

.. math::
    :label: second_vel

    \mathbf{u}^{k+1}\big|_w = \mathbf{u}^{**}\big|_w - \Delta{t}\cdot{c_k} \nabla{p}^{k+1}

At the walls, :math:`\mathbf{u}^*\big|_w=0` (no-slip condition for the first intermediate velocity) which when substituted into Eq. :eq:`first_vel` gives:

.. math::
    :label: second_vel2

    \mathbf{u}^{**}\big|_w=\Delta{t}\cdot{c_k}\nabla{p}^K

Therefore, the velocity at the wall in the current time-step is:

.. math::
    :label: wall_vel1

    \mathbf{u}^{k+1}\big|_w =\Delta{t}\cdot{c_k}\left(\nabla{p}^{k}-\nabla{p}^{k+1}\right)

Since for small time steps :math:`\nabla{p}^{k+1}\approx\nabla{p}^{k}` this results in:


.. math::
    :label: wall_vel2

    \mathbf{u}^{k+1}\big|_w\approx{0}

which ensures that the no-slip boundary condition is satisfied.

Pressure treatment
~~~~~~~~~~~~~~~~~~

The incompressibility condition :eq:`div-free` can be verified at the end of each sub-time step :math:`\nabla\mathbf{u}\cdot\mathbf{u}^{k+1} =0`  through the solving of a Poisson equation:

.. math::
    :label: poisson

    \nabla\cdot\nabla \tilde{p}^{k+1} = \frac{\nabla\cdot\mathbf{u}^{**}}{c_k\Delta{t}}

that provides the estimation of :math:`\tilde{p}^{k+1}` required to perform the pressure correction.


.. _spatial_discretisation:

Spatial discretisation
----------------------

Convective and viscous terms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Assuming that we have a uniform distribution of :math:`n_x` nodes :math:`x_i` on the domain :math:`[0, L_x]` with :math:`x_i=(i-1)\Delta{x}` for :math:`1\le{i}\le{n_x}`, the first derivative :math:`f'(x)` of the function :math:`f(x)` can be approximated by a finite difference scheme of the form:

.. math::
    :label: first-derivative

    \alpha{f}'_{i-1} + f'_i + \alpha{f}'_{i+1} = a\frac{f_{i+1}-f_{i-1}}{2\Delta{x}} + b\frac{f_{i+2}-f_{i-2}}{4\Delta{x}}

By choosing :math:`\alpha=1/3`, :math:`a=14/9`, :math:`b=1/9` this approximation is sixth-order accurate while having a so-called "quasi-spectral behaviour" due to its capabilities to represent accurately a wide range of scales. The compromise of the sixth-order accuracy has been chosen to maintain a compact formulation via the use of a Hermitian structure of the scheme with :math:`\alpha\ne{0}`. Even though this scheme is twice as expensive as a second-order scheme, in order to get the same solution with second-order scheme will require four to five times more mesh nodes.

Pressure
~~~~~~~~

Convective and diffusive terms are discretised using scheme :eq:`first-derivative` on a collocated mesh whereas a partially staggered mesh is used for the pressure treatment. To evaluate :math:`f'_{i+1/2}` of the first derivative at the staggered nodes by a half-mesh :math:`\Delta{x}/2`, the sixth-order finite-difference scheme can be expressed as:

.. math::
    :label: staggered-first-derivative

    \alpha{f}'_{i-1/2} + f'_{i+1/2} + \alpha{f}'_{i+3/2} = a\frac{f_{i+1}-f_{i-1}}{\Delta{x}} + b\frac{f_{i+2}-f_{i-1}}{3\Delta{x}}

with :math:`\alpha=9/62`, :math:`a=63/62` and :math:`b=17/62`. The spectral behaviour of this scheme is better than its collocated counterpart :eq:`first-derivative`.

.. tikz:: Arrangement of variables in 2D for a partially staggered grid.

    \draw[step=0.6cm,gray,very thin] (-0.5,-0.5) grid (2.9,2.9);    % Adjusted grid size
    \draw[step=1.2cm,gray,very thin] (-0.5,-0.5) grid (2.9,2.9);     % Adjusted grid size
    \foreach \x in {0,1.2,2.4}   
        \foreach \y in {0,1.2,2.4}   
            \fill (\x,\y) circle (2pt);
    \foreach \x in {0.6,1.8}     
        \foreach \y in {0.6,1.8} {    
            \fill[white] (\x,\y) circle (2pt); 
            \draw (\x,\y) circle (2pt);       
        }

    \draw[<->] (0,-0.5) -- (1.2,-0.5) node[midway,below] {\fontsize{2}{2}\selectfont $\Delta x$};   
    \draw[<->] (2.9,0) -- (2.9,1.2) node[midway,right] {\fontsize{2}{2}\selectfont $\Delta y$}; 

    \node[left] at (-0.5,0) {\fontsize{2}{2}\selectfont $j-1$};     
    \node[left] at (-0.5,0.6) {\fontsize{2}{2}\selectfont $j-\frac{1}{2}$};
    \node[left] at (-0.5,1.2) {\fontsize{2}{2}\selectfont $j$};
    \node[left] at (-0.5,1.8) {\fontsize{2}{2}\selectfont $j+\frac{1}{2}$};
    \node[left] at (-0.5,2.4) {\fontsize{2}{2}\selectfont $j+1$};

    \node[above] at (0,2.9) {\fontsize{2}{2}\selectfont $i-1$};    
    \node[above] at (0.6,2.9) {\fontsize{2}{2}\selectfont $i-\frac{1}{2}$};
    \node[above] at (1.2,2.9) {\fontsize{2}{2}\selectfont $i$};
    \node[above] at (1.8,2.9) {\fontsize{2}{2}\selectfont $i+\frac{1}{2}$};
    \node[above] at (2.4,2.9) {\fontsize{2}{2}\selectfont $i+1$};

    % Key
    \node[anchor=west] at (-1,-1.0) {\fontsize{2}{2}\selectfont $\bullet$ u, v \quad  $\circ$ p};


Assuming that :math:`f` is periodic over the domain :math:`[0,L_x]`, the discrete Fourier transform of the function :math:`f` can be expressed as:

.. math::
    :label: dft

    \hat{f}_l = \frac{1}{n_x}\sum_{i=1}^{n_x} f_i e^{-{i}k_xx_i} \quad \mathrm{for} \quad -n_x/2 \le l \le n_x/2-1

where :math:`l=\sqrt{-1}` and :math:`k_x=2\pi{l}/L_x` is the wave number. The inverse discrete Fourier transform is given by:

.. math::
    :label: idft

    f_i = \sum_{l=-n_x/2}^{n_x/2-1} \hat{f}_l e^{i{k_x}x_i}.

It can be shown that the Fourier coefficients :math:`\hat{f}'_l` associated with the approximation :eq:`first-derivative` are linked to the Fourier coefficients :math:`\hat{f}_l` given by :eq:`dft` by the simple spectral relation:

.. math::
    :label: spectral-relation

    \hat{f}'_l = l{k'_x}\hat{f}_l

where :math:`k'_x` is the modified wave number related to the actual wave number :math:`k_x` by

.. math::
    :label: modified-wave-number

    k'_x\Delta{x} = \frac{a\sin(k_x\Delta{x}) + (b/2)\sin(2k_x\Delta{x})}{1+2\alpha\cos(k_x\Delta{x})}

The concept of the modified wave number still holds in the staggered formulation, and the expression of :math:`k'_x` associated with the scheme :eq:`staggered-first-derivative` is given by:

.. math::
    :label: modified-wave-number-staggered

    k'_x\Delta{x} = \frac{2a\sin(k_x\Delta{x}/2) + (2b/3)\sin(3k_x\Delta{x}/2)}{1+2\alpha\cos(k_x\Delta{x})}

The well known principle of equivalence between multiplication in Fourier space and derivation/interpolation in the physical space is recalled here. This equivalence is exact, hence, the computation of a derivative in physical space using :eq:`staggered-first-derivative` with relevant boundary conditions must lead to the same result obtained with the use of :eq:`modified-wave-number-staggered` in spectral space.

.. _solving_the_poisson_equation:

Solving the Poisson equation
----------------------------

There are several numerical algorithms for solving Poisson's equations, which can be broadly classified into two categories: iterative solvers and direct solvers. x3d2 currently uses direct methods with iterative solvers planned in future versions. Among the direct methods, Fast Fourier Transform (FFT) based solvers are the most efficient.

For simplicity, a generic 3D Fourier transform can be defined as:

.. math::
    :label: 3d-dft

    \hat{p}_{lmn} = \frac{1}{n_xn_yn_z}\sum_{i}\sum_{j}\sum_{k} p_{ijk} W_x(k_xx_i)W_y(k_yy_j)W_z(k_zz_k)

with its inverse expression

.. math::
    :label: 3d-idft

    \hat{p}_{ijk} = \sum_{l}\sum_{m}\sum_{n} \hat{p}_{lmn} W_x(-k_xx_i)W_y(k_yy_j)W_z(k_zz_k)

where the sums, the base functions :math:`(W_x, W_y, W_z)` and the wave numbers :math:`(k_x, k_y, k_z)` can correspond to standard FFT (for periodic boundary conditions) or cosine FFT (for free-slip or :math:`\mathbf{u}`-Dirichlet/:math:`p`-Neumann boundary conditions) in their collocated or staggered versions. 3D direct :eq:`3d-dft` and inverse :eq:`3d-idft` can be performed with any efficient FFT routines available in scientific Fortran or C libraries. The first stage in solving the Poisson equation :eq:`poisson` consists in the computation of its right-hand side. After performing the relevant Fourier transform :eq:`3d-dft` to :math:`D=\nabla\cdot\mathbf{u}^{**}`, the solving of the Poisson equation consists in a single division of each Fourier mode :math:`\hat{D}_{lmn}` by a factor :math:`F_{lmn}` with

.. math::
    :label: poisson-solve

    \hat{\tilde{p}}^{k+1}_{lmn} = \frac{\hat{D}_{lmn}}{F_{lmn}}

where the expression of this factor depends on the mesh configuration. For instance, in the case of a partially staggered approach, the factor :math:`F_{lmn}` must take the mid-point interpolation into account through the use of a transfer function with the following form:

.. math::
    :label: staggered-factor

    F_{lmn} = -[(k'_xT_yT_z)^2 + (k'_yT_xT_y)^2 + (k'_zT_xT_y)^2]c_k\Delta{t}

where :math:`T_x(k_x\Delta{x})` is the transfer function related to the wave number :math:`k_x` by

.. math::
    :label: transfer-function

    T_x(k_x\Delta{x}) = \frac{2a\cos(k_x\Delta{x}/2) + (2b/3)\cos(3k_x\Delta{x}/2)}{1+2\alpha\cos(k_x\Delta{x})}


Stretched mesh in one direction
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The pressure discretisation described so far is only valid for a regular mesh in three spatial directions. To overcome this difficulty a modification of the Poisson solver is proposed which is based on a specific function mapping and expressed using only few Fourier modes. This approach preserves the spectral and non-iterative nature of the pressure treatment without significant loss of accuracy. 

For simplicity consider a one-dimensional problem where :math:`y` is the physical coordinate and :math:`s` is the computational coordinate:

.. math::
    :label: mapping

    y = h(s), \quad 0\le{s}\le{1}, 0\le{y}\le{L_y}

where :math:`h(s)` is the mapping from equally spaced coordinate :math:`s` to the stretched physical coordinate :math:`y`. The derivatives with respect to :math:`y` can be estimated using the chain rule, where the first derivative is given by:

.. math::
    :label: first-derivative-mapping

    \frac{\partial{f}}{\partial{y}} = \frac{\partial{f}}{\partial{s}}\frac{\partial{s}}{\partial{y}} = \frac{1}{h'(s)}\frac{\partial{f}}{\partial{s}}

and the second derivative is given by:

.. math::
    :label: second-derivative-mapping

    \frac{\partial^2{f}}{\partial{y}^2} = \frac{\partial^2{f}}{\partial{s^2}}\left(\frac{\partial{s}}{\partial{y}}\right)^2 + \frac{\partial{f}}{\partial{s}}\frac{\partial^2{s}}{\partial{y^2}} = \frac{1}{h'(s)^2}\frac{\partial^2{f}}{\partial{s}^2} - \frac{h''(s)}{h'(s)^3}\frac{\partial{f}}{\partial{s}}

Expressed in physical space, these rules can be used to implement schemes like :eq:`first-derivative` and :eq:`staggered-first-derivative` where the finite differences are performed on the regular coordinate :math:`s` (instead of :math:`x`). 

The main difficulty is in the treatment of the Poisson equation that requires similar operations in the spectral space. Here the metric :math:`1/h'` is expressed with only three Fourier modes in spectral space:

.. math::
    :label: truncated-metric

    \frac{1}{h'} = \frac{1}{L_y}\left\{ \frac{\alpha}{\pi} + \frac{1}{\pi\beta}\sin^2(\pi(\gamma{s} + \delta))\right\}  = \frac{1}{L_y}\left\{ \frac{\alpha}{\pi} + \frac{1}{2\pi\beta}\left[1-\frac{e^{l2\pi(\gamma{s}+\delta)} + e^{-l2\pi(\gamma{s}+\delta)}}{2}\right]\right\}

so that the mapping :eq:`mapping` can be written as:

.. math::
    :label: fourier-mapping-metric

    \begin{align}
    h &= \frac{L_y\sqrt{\beta}}{\gamma\sqrt\alpha\sqrt{\alpha\beta+1}}\left\{\tan^{-1}\left[\frac{\sqrt{\alpha\beta+1}\tan(\pi(\gamma{s}+\delta))}{\sqrt\alpha\sqrt\beta}\right] \right. \\
    &+ \left. \pi\left[H\left(s-\frac{1-2\delta}{2\gamma}\right) + H\left(s-\frac{3-2\delta}{2\gamma}\right)\right] -\tan^{-1}\left[\frac{\sqrt{\alpha\beta+1} +\tan(\pi\delta)}{\sqrt\alpha\sqrt\beta}\right]   \right\}
    \end{align}

where :math:`H` is the Heaviside step function. This mapping preserves the accuracy while avoiding expensive computation of a full convolution and ensuring the strict physical/spectral equivalence.

* :math:`\alpha=0`, :math:`\gamma=1` and :math:`\delta=0` the mapping leads to refinement in the centre of an infinite domain
* :math:`\alpha\ne{0}`, :math:`\gamma=1` and :math:`\delta=0` leads to refinement in the centre of a finite domain
* :math:`\gamma=1` and :math:`\delta=1/2` leads to refinement near the boundaries for a finite domain (not compatible with periodic boundary conditions) because :math:`1/h'` is not periodic over :math:`L_y`
* :math:`\gamma=1/2` and :math:`\delta=1/2` leads to refinement near the bottom boundary only for a finite domain

It can be deduced that the three coefficients of the metric :eq:`fourier-mapping-metric` are non-zero with

.. math::
    :label: metric-coefficients

    \alpha = \frac{1}{L_y}\left(\frac{\alpha}{\pi} + \frac{1}{2\pi\beta}\right), \quad \hat{a}_1 = \hat{a}_{-1}  = -\frac{1}{L_y}\left(\frac{\cos{2\pi\delta}}{4\pi\beta}\right)

for :math:`\gamma=1` and :math:`\delta=0` or :math:`1/2`. The main advantage of this compact expression in spectral space is that the convolution of the metric by the first derivation with respect to the regular coordinate :math:`s` requires only :math:`3n_y` multiplications.

To solve the Poisson equation :eq:`poisson` (using 3D Fourier transforms :eq:`3d-dft` and :eq:`3d-idft` where :math:`y` needs to be substituted by :math:`s` for the :math:`y`-stretched approach) the counterpart of the integration scheme :eq:`poisson-solve` becomes


.. math::
    :label: poisson-solve-stretched

    \hat{\tilde{\mathbf{p}}}_{ln}^{k+1} = \mathbf{B}^{-1}\widehat{\mathbf{D}}_{ln}

where :math:`\hat{\tilde{\mathbf{p}}}_{ln}^{k+1}` and :math:`\widehat{\mathbf{D}}_{ln}` are :math:`n_y` vectors of components of :math:`\hat{\tilde{p}}_{ln}^{k+1}` and :math:`\widehat{D}_{ln}` and :math:`\mathbf{B}` is a :math:`n_y \times n_y` pentadiagonal matrix of components. For the partially staggered case these components are:

.. math::
    :label: pentadiagonal-matrix

    &b_{m,m-2} = -\hat{a}_1^2T_x^2T_z^2k'_{m-1}k'_{m-2} \\
    &b_{m,m-1} = -\hat{a}_0\hat{a}_1^2T_x^2T_z^2k'_{m-1}(k'_m + k'_{m-1}) \\
    &b_{m,m} = -(k'_xT_yT_z)^2 - (k'_zT_yT_z)^2 -\hat{a}^2_0T_x^2T_z^2{k'_m}^2 - \hat{a}_1\hat{a}_{-1}T_x^2T_z^2k'_m(k'_{m+1} + k'_{m-1}) \\
    &b_{m,m+1} = -\hat{a}_0\hat{a}_1T_x^2T_z^2k'_{m+1}(k'_m + k'_{m+1}) \\
    &b_{m,m+2} = -\hat{a}_{-1}^2T_x^2T_z^2k'_{m+1}k'_{m+2}

where the :math:`k'_m` are the modified wave numbers from relation like :eq:`modified-wave-number` or :eq:`modified-wave-number-staggered` based on the computational coordinate :math:`s` instead of :math:`x`. 

The above matrix is diagonal for a regular :math:`y`-coordinate (with :math:`a_1=a_{-1}=0`) so that the simplified expression :eq:`transfer-function` can be recovered. In the other cases, the computation of pressure nodes :math:`\hat{\tilde{\mathbf{p}}}_{ln}^{k+1}` requires inverting :math:`n_x \times n_y` linear systems based on :math:`n_y\times{n_y}` pentadiagonal matrices. The corresponding computational cost is proportional to :math:`n_x\times{n_y}\times{n_z}` so that the solver Poisson can be direct without any iterative process.

In terms of computational cost, solving the Poisson equation directly requires both a forward and an inverse 3D FFT. For a completely regular mesh in three spatial dimensions, these two FFT operations constitute the majority of the computational expense for the Poisson stage, accounting for about 10% of the total computational effort required to solve the Navier-Stokes equations. When dealing with meshes that have one stretched direction, the cost of ensuring incompressibility increases but still represents about 15% of the overall computational cost for a given simulation.

Although using Fourier transforms for pressure is highly suitable for periodic or free-slip boundary conditions, it is less ideal for no-slip or open boundary conditions. In these cases the pressure must be expressed using cosine Fourier transforms, assuming that homogeneous Neumann conditions are met. This assumption introduces an error that is only second-order accurate in space.

.. _tridiagonal_systems:

Tridiagonal systems
-------------------

A tridiagonal system is a linear system of equations where the non-zero coefficients are concentrated along the main diagonal, the subdiagonal, and the superdiagonal. These systems frequently occur in high-order compact finite difference schemes and take the form :math:`Ax=d`:

.. math::
    :label: tdsops1

    a_iu_{i-1} + b_iu_i + c_iu_{i+1} = d_i \quad i=0,1,\dots,N-1,

where :math:`a_0=c_{N-1}=0`.

Tridiagonal systems solver algorithms
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The Thomas algorithm is a well-known method for solving tridiagonal systems of equations. It is a specialised form of Gaussian elimination that involves a forward pass to eliminate the lower diagonal elements :math:`a_i` of the tridiagonal matrix by adding a multiple of the row above, followed by a backward substitution pass using the modified upper diagonal coefficients :math:`c_i`. Even though the Thomas algorithm is inherently serial, when solving a batch of tridiagonal systems it is possible to parallelise the overall operation by assigning many tridiagonal systems per MPI rank. 


However, the Thomas algorithm requires the entire length of a tridiagonal system to be contained within a single rank. Consequently, this strategy necessitates MPI all-to-all communication to rearrange the domain decomposition, ensuring that batches of tridiagonal systems are oriented in the :math:`x-`, :math:`y-`, and :math:`z-` directions. This ensures that the entire length of each tridiagonal system is always within a single rank, while individual systems are distributed across all available ranks. As a result, the derivations and interpolations can be performed in each direction using the Thomas algorithm.

The fundamental difference between a serial algorithms like Thomas algorithm and distributed-memory algorithms is that the distributed algorithms divide the individual systems into multiple subdomains. This enables localised communication (e.g. neighbour-to-neighbour exchanges) instead of global dependencies.  For example, a 3D decomposition strategy can distribute subdomains across ranks while maintaining static decomposition states, avoiding frequent reconfiguration. Crucially, these approaches eliminate the need for MPI all-to-all communications across the entire domain.

Many tridiagonal algorithms face performance challenges, particularly along the :math:`x` direction, due to inefficient memory access patterns in Cartesian data structures. There are two main characteristics in tridiagonal matrix algorithms that are crucial for computational performance:

1. The algorithms are heavily bandwidth-bound because the FLOP requirements of these algorithms are very low compared to the data movement requirements.
2. The backward and forward passes of the algorithm require accessing the same data twice in a short time interval. This can be optimised by using CPU cache or GPU shared memory to store intermediate states between the forward and backward passes.

x3d2's solution: DistD2-TDS
~~~~~~~~~~~~~~~~~~~~~~~~~~~

To address the challenges of solving tridiagonal systems efficiently in distributed-memory environments on both CPUs and GPUs, x3d2 uses a novel algorithm called `DistD2-TDS` (see :cite:`akkurt_cpc_24` for more details). This algorithm is based on a specialised data structure that:

1. Improves data locality and minimises data movements via cache blocking and kernel fusion strategies.
2. Enables a contiguous data access pattern, resulting in efficient utilisation of the available memory bandwidth.
3. Supports vectorisation on CPUs and thread-level parallelisation on GPUs for improved performance.

Modern CPUs have large enough caches to store the entire tridiagonal matrix, which can be accessed quickly. However, GPU caches are not large enough, and kernel fusion is the only option to reduce data movement on GPUs. To enable a linear and predictive memory access pattern regardless of the spatial direction of the tridiagonal systems, DistD2-TDS data structure subdivides the computational domain into groups of individual tridiagonal systems. These groups are tightly packed in memory to ensure data continuity. This arrangement enables vectorisation on CPUs and thread-level parallelism on GPUs, as the :math:`n^{th}` entries of all tridiagonal systems within a group are stored next to each other in memory. Therefore, the sequential operations in the algorithms, as we apply the forward and backward passes, can be concurrently executed for :math:`SZ` systems at once per core on a CPU or per SM on a GPU.

Example: In a domain of size :math:`32 \times 8 \times 4` shown below; a single group consist of 4 individual tridiagonal systems, resulting in 8 groups in total. For a Cartesian mesh with :math:`n_x`, :math:`n_y`, and :math:`n_z` the data is arranged as  :math:`SZ, n_x, n_y \cdot n_z / SZ` for an :math:`x`-directional layout where :math:`SZ` is the number of tridiagonal systems in a group.

.. tikz::  x3d2 data structure for an :math:`x`-directional tridiagonal system. Data continuity in memory is in column-major order.

   % first diagram
   \begin{scope}[rotate around y=90]

   % axes
   \draw[thick,->] (0,0,0) -- (10,0,0) node[anchor=north east, xshift=10]{$x$}; % x-axis
   \draw[thick,->] (0,0,0) -- (0,3,0) node[anchor=south]{$y$}; % y-axis
   \draw[thick,->] (0,0,0) -- (0,0,4) node[anchor=east, xshift=12]{$z$}; % z-axis

   % parameters
   \def\dx{0.2}
   \def\dy{0.1}
   \def\dz{0.1}
   \def\nx{40} 
   \def\ny{7}
   \def\nz{3}

   % front face
   \foreach \x in {0,...,40} {
     \foreach \y in {0,...,3} {
       \pgfmathsetmacro\k{\x*2.5}
       \fill [color=red!\k!blue] (\x*\dx, \y*\dy, 0) circle (1pt);
     }
   }

   % boundary points
   \foreach \x in {0,\nx} {
     \foreach \y in {0,...,7} {
       \foreach \z in {0,...,3} {
         \pgfmathsetmacro\k{\x*2.5}
         \fill [color=red!\k!blue] (\x*\dx, \y*\dy, \z*\dz) circle (1pt);
       }
     }
   }

   % box structure
   \path (0,0,0) coordinate (A)
         (0,0,\nz*\dz) coordinate (B)
         (0,\ny*\dy,0) coordinate (C)
         (0,\ny*\dy,\nz*\dz) coordinate (D)
         (\nx*\dx,0,0) coordinate (E)
         (\nx*\dx,0,\nz*\dz) coordinate (F)
         (\nx*\dx,\ny*\dy,0) coordinate (G)
         (\nx*\dx,\ny*\dy,\nz*\dz) coordinate (H);

   % draw the box
   \draw (A)--(B)--(D)--(C)--(A);
   \draw (E)--(F)--(H)--(G)--(E);
   \draw (A)--(E);
   \draw (B)--(F);
   \draw (C)--(G);
   \draw (D)--(H);

   %% internal lines
   %\foreach \i in {1,2,3} {
   %  \draw (0, \i*\dy, 0) -- (\nx*\dx, \i*\dy, 0);
   %}
   \end{scope}

   % second diagram (shifted to the right)
   \begin{scope}[xshift=5cm]  % Adjust spacing between the two diagrams

   % parameters
   \def\dx{0.1}
   \def\dy{0.1}
   \def\dz{0.1}
   \def\nx{31}
   \def\ny{7}
   \def\nz{3}

   % colours
   \definecolor{redw}{rgb}{1,0.45,0.45}
   \definecolor{bluew}{rgb}{0.45,0.45,1}

   % sub-groups
   \foreach \i in {1,...,8} {
     \foreach \x in {0,...,31} {
       % group labels
       \node[](Group\i) at (-0.3, 1.5*8*4*\dy - 0.2 - \i*4*\dy - \i*0.1) {\i};

       % grid points
       \foreach \y in {0,...,3} {
         \foreach \z in {0} {
           \pgfmathsetmacro\k{\x*2.5}
           \draw [color=red!\k!blue, mark=*, mark size=1]
                 plot coordinates {(\x*\dx, \y*\dy + \i*4*\dy + \i*0.1, \z*\dz)};
         }
       }
     }
   }
   \end{scope}

Typically, for a double-precision simulation on CPUs, :math:`SZ=8` as vector registers (512 bits) handle 8 double-precision FLOPs per cycle. For GPUs, :math:`SZ=32` as a single streaming multiprocessor (SM) typically has 32 double-precision cores in total, where each core is effectively assigned an individual tridiagonal system.

Tridiagonal systems resulting from high-order compact finite-difference schemes are always diagonally dominant, and the algorithm used in x3d2 takes advantage of this characteristic to minimise communication requirements between ranks. This involves dividing a batch of tridiagonal systems into multiple domains, where the subdomains are located across multiple ranks in a distributed memory environment. 

Apart from the specialised data structure in DistD2, another novel aspect of this algorithm is that it fuses the RHS construction based on discretised equations such as Eq. :eq:`first-derivative` with the forward pass in the decoupling phase of the algorithm -- this minimises data movement requirements.

After the decoupling phase of the algorithm, the next step is to solve the reduced system. By exploiting diagonal dominance, certain entries in the reduced system  can be eliminated without any loss in numerical accuracy, leading to independent :math:`2 \times 2` systems that are coupled across MPI boundaries. This step is the main deviation from other tridiagonal solvers that use a PCR-type strategy to obtain the solution of the reduced system on a single rank. In the DistD2 algorithm, the reduced systems are not coupled across the entire domain, but only across 2 individual ranks across an MPI boundary. These smaller systems are then solved efficiently using local communication between neighbouring ranks. This single-step MPI communication involving two neighbouring ranks is a significant advantage of the DistD2 algorithm over existing strategies. Finally, the substitution phase of the algorithm requires a simple algebraic substitution and does not require any further MPI communications.

2D Domain Decomposition for FFT
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The current implementation of x3d2 uses a direct approach, namely Fast Fourier Transform (FFT), to solve the Poisson equation. Spectral codes often involve performing FFTs along global mesh lines. There are two approaches to performing such computations on distributed-memory systems: 

1. Develop distributed algorithms (such as a parallel tridiagonal solver or parallel FFT working on distributed data).
2. Dynamically redistribute (transpose) data among processors to apply serial algorithms in local memory.

The second approach is often preferred due to its simplicity. Many applications have implemented this idea using 1D domain decomposition. However, 1D decomposition has some limitations. For example, for a cubic mesh size of :math:`N^3`, the maximum number of processors :math:`N_{\mathrm{proc}}` that can be used in a 1D decomposition is :math:`N`, as each slab must contain at least one plane. For a cubic mesh of 1 billion points, this constraint is :math:`N_{\mathrm{proc}} \lt 10^4`. This limitation can be overcome by using 2D domain decomposition. While a 1D decomposition algorithm swaps between two states, a 2D decomposition requires traversing three different states using four global transpositions to complete a cycle. The swapping between states can be achieved using the `MPI_ALLTOALL(V)` library.

.. tikz:: 1D domain decomposition (left), 2D domain decomposition (middle), 3D domain decomposition (right).
   :align: center
   :libs: arrows, shapes.geometric, shadows, fit, patterns, automata, quotes, arrows.meta, decorations.pathreplacing, calligraphy

   \newif\ifcuboidshade
   \newif\ifcuboidemphedge
 
   \tikzset{
     cuboid/.is family,
     cuboid,
     shiftx/.initial=0,
     shifty/.initial=0,
     dimx/.initial=3,
     dimy/.initial=3,
     dimz/.initial=3,
     scale/.initial=1,
     densityx/.initial=1,
     densityy/.initial=1,
     densityz/.initial=1,
     rotation/.initial=0,
     anglex/.initial=0,
     angley/.initial=90,
     anglez/.initial=225,
     scalex/.initial=1,
     scaley/.initial=1,
     scalez/.initial=0.5,
     front/.style={draw=black,fill=white},
     top/.style={draw=black,fill=white},
     right/.style={draw=black,fill=white},
     shade/.is if=cuboidshade,
     shadecolordark/.initial=black,
     shadecolorlight/.initial=white,
     shadeopacity/.initial=0.15,
     shadesamples/.initial=16,
     emphedge/.is if=cuboidemphedge,
     emphstyle/.style={thick},
   }
 
   \newcommand{\tikzcuboidkey}[1]{\pgfkeysvalueof{/tikz/cuboid/#1}}
 
   % commands
   \newcommand{\tikzcuboid}[1]{
     \tikzset{cuboid,#1} % Process Keys passed to command
     \pgfmathsetlengthmacro{\vectorxx}{\tikzcuboidkey{scalex}*cos(\tikzcuboidkey{anglex})*28.452756}
     \pgfmathsetlengthmacro{\vectorxy}{\tikzcuboidkey{scalex}*sin(\tikzcuboidkey{anglex})*28.452756}
     \pgfmathsetlengthmacro{\vectoryx}{\tikzcuboidkey{scaley}*cos(\tikzcuboidkey{angley})*28.452756}
     \pgfmathsetlengthmacro{\vectoryy}{\tikzcuboidkey{scaley}*sin(\tikzcuboidkey{angley})*28.452756}
     \pgfmathsetlengthmacro{\vectorzx}{\tikzcuboidkey{scalez}*cos(\tikzcuboidkey{anglez})*28.452756}
     \pgfmathsetlengthmacro{\vectorzy}{\tikzcuboidkey{scalez}*sin(\tikzcuboidkey{anglez})*28.452756}

     \begin{scope}[xshift=\tikzcuboidkey{shiftx}, yshift=\tikzcuboidkey{shifty}, scale=\tikzcuboidkey{scale},
                   rotate=\tikzcuboidkey{rotation}, x={(\vectorxx,\vectorxy)}, y={(\vectoryx,\vectoryy)}, z={(\vectorzx,\vectorzy)}]

        \pgfmathsetmacro{\steppingx}{1/\tikzcuboidkey{densityx}}
        \pgfmathsetmacro{\steppingy}{1/\tikzcuboidkey{densityy}}
        \pgfmathsetmacro{\steppingz}{1/\tikzcuboidkey{densityz}}

        \newcommand{\dimx}{\tikzcuboidkey{dimx}}
        \newcommand{\dimy}{\tikzcuboidkey{dimy}}
        \newcommand{\dimz}{\tikzcuboidkey{dimz}}

        \pgfmathsetmacro{\secondx}{2*\steppingx}
        \pgfmathsetmacro{\secondy}{2*\steppingy}
        \pgfmathsetmacro{\secondz}{2*\steppingz}

        \foreach \x in {\steppingx,\secondx,...,\dimx}
        { \foreach \y in {\steppingy,\secondy,...,\dimy}
          { \pgfmathsetmacro{\lowx}{(\x-\steppingx)}
            \pgfmathsetmacro{\lowy}{(\y-\steppingy)}
            \filldraw[cuboid/front] (\lowx,\lowy,\dimz) -- (\lowx,\y,\dimz) -- (\x,\y,\dimz) -- (\x,\lowy,\dimz) -- cycle;
          }
        }
        \foreach \x in {\steppingx,\secondx,...,\dimx}
        { \foreach \z in {\steppingz,\secondz,...,\dimz}
          { \pgfmathsetmacro{\lowx}{(\x-\steppingx)}
            \pgfmathsetmacro{\lowz}{(\z-\steppingz)}
            \filldraw[cuboid/top] (\lowx,\dimy,\lowz) -- (\lowx,\dimy,\z) -- (\x,\dimy,\z) -- (\x,\dimy,\lowz) -- cycle;
              }
        }
          \foreach \y in {\steppingy,\secondy,...,\dimy}
        { \foreach \z in {\steppingz,\secondz,...,\dimz}
          { \pgfmathsetmacro{\lowy}{(\y-\steppingy)}
            \pgfmathsetmacro{\lowz}{(\z-\steppingz)}
            \filldraw[cuboid/right] (\dimx,\lowy,\lowz) -- (\dimx,\lowy,\z) -- (\dimx,\y,\z) -- (\dimx,\y,\lowz) -- cycle;
          }
        }
        \ifcuboidemphedge
          \draw[cuboid/emphstyle] (0,\dimy,0) -- (\dimx,\dimy,0) -- (\dimx,\dimy,\dimz) -- (0,\dimy,\dimz) -- cycle;%
          \draw[cuboid/emphstyle] (0,\dimy,\dimz) -- (0,0,\dimz) -- (\dimx,0,\dimz) -- (\dimx,\dimy,\dimz);%
          \draw[cuboid/emphstyle] (\dimx,\dimy,0) -- (\dimx,0,0) -- (\dimx,0,\dimz);%
        \fi
 
        \ifcuboidshade
          \pgfmathsetmacro{\cstepx}{\dimx/\tikzcuboidkey{shadesamples}}
          \pgfmathsetmacro{\cstepy}{\dimy/\tikzcuboidkey{shadesamples}}
          \pgfmathsetmacro{\cstepz}{\dimz/\tikzcuboidkey{shadesamples}}
          \foreach \s in {1,...,\tikzcuboidkey{shadesamples}}

          { \pgfmathsetmacro{\lows}{\s-1}
            \pgfmathsetmacro{\cpercent}{(\lows)/(\tikzcuboidkey{shadesamples}-1)*100}
            \fill[opacity=\tikzcuboidkey{shadeopacity},color=\tikzcuboidkey{shadecolorlight}!\cpercent!\tikzcuboidkey{shadecolordark}] (0,\s*\cstepy,\dimz) -- (\s*\cstepx,\s*\cstepy,\dimz) -- (\s*\cstepx,0,\dimz) -- (\lows*\cstepx,0,\dimz) -- (\lows*\cstepx,\lows*\cstepy,\dimz) -- (0,\lows*\cstepy,\dimz) -- cycle;
            \fill[opacity=\tikzcuboidkey{shadeopacity},color=\tikzcuboidkey{shadecolorlight}!\cpercent!\tikzcuboidkey{shadecolordark}] (0,\dimy,\s*\cstepz) -- (\s*\cstepx,\dimy,\s*\cstepz) -- (\s*\cstepx,\dimy,0) -- (\lows*\cstepx,\dimy,0) -- (\lows*\cstepx,\dimy,\lows*\cstepz) -- (0,\dimy,\lows*\cstepz) -- cycle;
            \fill[opacity=\tikzcuboidkey{shadeopacity},color=\tikzcuboidkey{shadecolorlight}!\cpercent!\tikzcuboidkey{shadecolordark}] (\dimx,0,\s*\cstepz) -- (\dimx,\s*\cstepy,\s*\cstepz) -- (\dimx,\s*\cstepy,0) -- (\dimx,\lows*\cstepy,0) -- (\dimx,\lows*\cstepy,\lows*\cstepz) -- (\dimx,0,\lows*\cstepz) -- cycle;
          }
        \fi
     \end{scope}
   }
 
   \makeatother
 
   % simple cross
   \tikzset{cross/.style={cross out, draw=black, minimum size=2*(#1-\pgflinewidth), inner sep=0pt, outer sep=0pt},
   % default radius will be 1pt.
   cross/.default={1pt}}
 
   % 1D
   \tikzcuboid{%
    shiftx=0cm,%
    shifty=0cm,%
    scale=1,%
    rotation=0,%
    densityx=2,%
    densityy=2,%
    densityz=2,%
    dimx=5,%
    dimy=1,%
    dimz=5,%
    scalex=0.95,%
    scaley=1,%
    scalez=1,%
    anglex=0,%
    angley=90,%
    anglez=225,%
    front/.style={draw=red!50!black,fill=red!50!white},%
    top/.style={draw=red!50!black,fill=red!50!white},%
    right/.style={draw=red!50!black,fill=red!50!white},%
    emphedge=false,%
   }

   \begin{scope}[shift={(0,1.5,0)}]
       \tikzcuboid{%
       front/.style={draw=green!50!black,fill=green!50!white},%
       top/.style={draw=green!50!black,fill=green!50!white},%
       right/.style={draw=green!50!black,fill=green!50!white},%
       };
   \end{scope}

   \begin{scope}[shift={(0,3,0)}]
       \tikzcuboid{%
       front/.style={draw=blue!50!black,fill=blue!50!white},%
       top/.style={draw=blue!50!black,fill=blue!50!white},%
       right/.style={draw=blue!50!black,fill=blue!50!white},%
       }
   \end{scope}

   \begin{scope}[shift={(0,4.5,0)}]
       \tikzcuboid{%
       front/.style={draw=yellow!50!black,fill=yellow!50!white},%
       top/.style={draw=yellow!50!black,fill=yellow!50!white},%
       right/.style={draw=yellow!50!black,fill=yellow!50!white},%
       }
   \end{scope}

   \draw[thick,->] (-2.5,4.5,0) -- (-1.5,4.5,0) node[anchor=north east, font=\Large]{$x$};
   \draw[thick,->] (-2.5,4.5,0) -- (-2.5,5.5,0) node[anchor=north west, font=\Large]{$y$};
   \draw[thick,->] (-2.5,4.5,0) -- (-2.5,4.5,1) node[anchor=east, font=\Large]{$z$};

   % 2D
   \begin{scope}[xshift=10cm]
     \tikzcuboid{%
      shiftx=0cm,%
      shifty=0cm,%
      scale=1,%
      rotation=0,%
      densityx=2,%
      densityy=2,%
      densityz=2,%
      dimx=1,%
      dimy=1,%
      dimz=5,%
      scalex=1,%
      scaley=1,%
      scalez=0.8,%
      anglex=0,%
      angley=90,%
      anglez=225,%
      front/.style={draw=red!50!black,fill=red!50!white},%
      top/.style={draw=red!50!black,fill=red!50!white},%
      right/.style={draw=red!50!black,fill=red!50!white},%
      emphedge=false,%
     }

     \begin{scope}[shift={(0,1.5,0)}]
         \tikzcuboid{%
         front/.style={draw=green!50!black,fill=green!50!white},%
         top/.style={draw=green!50!black,fill=green!50!white},%
         right/.style={draw=green!50!black,fill=green!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(0,3,0)}]
         \tikzcuboid{%
         front/.style={draw=blue!50!black,fill=blue!50!white},%
         top/.style={draw=blue!50!black,fill=blue!50!white},%
         right/.style={draw=blue!50!black,fill=blue!50!white},%
         }
     \end{scope}

     \begin{scope}[shift={(0,4.5,0)}]
         \tikzcuboid{%
         front/.style={draw=yellow!50!black,fill=yellow!50!white},%
         top/.style={draw=yellow!50!black,fill=yellow!50!white},%
         right/.style={draw=yellow!50!black,fill=yellow!50!white},%
         }
     \end{scope}
   
     \begin{scope}[shift={(1.5,0,0)}]
         \tikzcuboid{%
         front/.style={draw=lime!50!black,fill=lime!50!white},%
         top/.style={draw=lime!50!black,fill=lime!50!white},%
         right/.style={draw=lime!50!black,fill=lime!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(1.5,1.5,0)}]
         \tikzcuboid{%
         front/.style={draw=magenta!50!black,fill=magenta!50!white},%
         top/.style={draw=magenta!50!black,fill=magenta!50!white},%
         right/.style={draw=magenta!50!black,fill=magenta!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(1.5,3,0)}]
         \tikzcuboid{%
         front/.style={draw=olive!50!black,fill=olive!50!white},%
         top/.style={draw=olive!50!black,fill=olive!50!white},%
         right/.style={draw=olive!50!black,fill=olive!50!white},%
         }
     \end{scope}

     \begin{scope}[shift={(1.5,4.5,0)}]
         \tikzcuboid{%
         front/.style={draw=pink!50!black,fill=pink!50!white},%
         top/.style={draw=pink!50!black,fill=pink!50!white},%
         right/.style={draw=pink!50!black,fill=pink!50!white},%
         }  % #brown%teal%violet%lime%pink%purple
     \end{scope}
   
     \begin{scope}[shift={(3,0,0)}]
         \tikzcuboid{%
         front/.style={draw=white!50!black,fill=white!50!white},%
         top/.style={draw=white!50!black,fill=white!50!white},%
         right/.style={draw=white!50!black,fill=white!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(3,1.5,0)}]
         \tikzcuboid{%
         front/.style={draw=orange!50!black,fill=orange!50!white},%
         top/.style={draw=orange!50!black,fill=orange!50!white},%
         right/.style={draw=orange!50!black,fill=orange!50!white},%
         };    %purple
     \end{scope}

     \begin{scope}[shift={(3,3,0)}]
         \tikzcuboid{%
         front/.style={draw=gray!50!black,fill=gray!50!white},%
         top/.style={draw=gray!50!black,fill=gray!50!white},%
         right/.style={draw=gray!50!black,fill=gray!50!white},%
         }
     \end{scope}

     \begin{scope}[shift={(3,4.5,0)}]
         \tikzcuboid{%
         front/.style={draw=cyan!50!black,fill=cyan!50!white},%
         top/.style={draw=cyan!50!black,fill=cyan!50!white},%
         right/.style={draw=cyan!50!black,fill=cyan!50!white},%
         }
     \end{scope}

   \end{scope}

   % 3D
   \begin{scope}[xshift=20cm]
     \tikzcuboid{%
      shiftx=0cm,%
      shifty=-1cm,%
      scale=1,%
      rotation=0,%
      densityx=2,%
      densityy=2,%
      densityz=2,%
      dimx=3,%
      dimy=3,%
      dimz=1,%
      scalex=1,%
      scaley=1,%
      scalez=1,%
      anglex=0,%
      angley=90,%
      anglez=225,%
      front/.style={draw=gray!50!black,fill=gray!50!white},%
      top/.style={draw=gray!50!black,fill=gray!50!white},%
      right/.style={draw=gray!50!black,fill=gray!50!white},%
      emphedge=false,%
     }

     \begin{scope}[shift={(0,4,0)}]
         \tikzcuboid{%
         front/.style={draw=green!50!black,fill=green!50!white},%
         top/.style={draw=green!50!black,fill=green!50!white},%
         right/.style={draw=green!50!black,fill=green!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(4,0,0)}]
         \tikzcuboid{%
         front/.style={draw=red!50!black,fill=red!50!white},%
         top/.style={draw=red!50!black,fill=red!50!white},%
         right/.style={draw=red!50!black,fill=red!50!white},%
         }
     \end{scope}

     \begin{scope}[shift={(4,4,0)}]
         \tikzcuboid{%
         front/.style={draw=yellow!50!black,fill=yellow!50!white},%
         top/.style={draw=yellow!50!black,fill=yellow!50!white},%
         right/.style={draw=yellow!50!black,fill=yellow!50!white},%
         }
     \end{scope}
   
     \begin{scope}[shift={(0,0,3)}]
         \tikzcuboid{%
         front/.style={draw=cyan!50!black,fill=cyan!50!white},%
         top/.style={draw=cyan!50!black,fill=cyan!50!white},%
         right/.style={draw=cyan!50!black,fill=cyan!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(0,4,3)}]
         \tikzcuboid{%
         front/.style={draw=orange!50!black,fill=orange!50!white},%
         top/.style={draw=orange!50!black,fill=orange!50!white},%
         right/.style={draw=orange!50!black,fill=orange!50!white},%
         };
     \end{scope}

     \begin{scope}[shift={(4,0,3)}]
         \tikzcuboid{%
         front/.style={draw=lime!50!black,fill=lime!50!white},%
         top/.style={draw=lime!50!black,fill=lime!50!white},%
         right/.style={draw=lime!50!black,fill=lime!50!white},%
         }
     \end{scope}

     \begin{scope}[shift={(4,4,3)}]
         \tikzcuboid{%
         front/.style={draw=blue!50!black,fill=blue!50!white},%
         top/.style={draw=blue!50!black,fill=blue!50!white},%
         right/.style={draw=blue!50!black,fill=blue!50!white},%
         }
     \end{scope}

   \end{scope}
   
2D domain decomposition is widely used for spectral codes, particularly those compatible with implicit schemes in space. This method allows for efficient parallelization by dividing the computational domain into smaller subdomains, each handled by a separate processor. For a simulation with a cubic mesh of size :math:`N^3`, up to :math:`N^2` processors can be used, significantly increasing scalability compared to 1D decomposition.

x3d2 uses the `2DECOMP&FFT <https://2decomp-fft.github.io/>`_ library for 2D decomposition and FFT (see :cite:`li_cray_10` for more details) when using OpenMP as a backend (for NVIDIA GPUs it uses `cuFFT <https://developer.nvidia.com/cufft>`_). One of the key advantages of using the 2DECOMP&FFT library is that it does not require modifications to the existing derivative and interpolation subroutines, making it easier to implement. Additionally, this approach utilises customised global `MPI_ALLTOALL(V)` transpositions to redistribute data among processors. Although communication overhead can range from 30% to 80% of the total computational time, with up to 70 transpositions per time step, the overall efficiency and scalability of the simulations are greatly enhanced.

References
----------

.. bibliography::
    :cited:
    :style: unsrt