# Spectral Graph Wavelet Transform


> The spectral graph wavelet transform (SGWT) ([Hammond et al., 2011](https://www.sciencedirect.com/science/article/pii/S1063520310000552)) defines the wavelet transform for weighted graphs. The key point of the SGWT is to design a real-valued kernel function $g(\lambda)$ that acts as a bandpass filter in the Fourier domain. With the kernel functions $g(s\lambda)$ at discrete scales $s$, we are able to create wavelets that are not only localized in the graph spectrum, but also localized in the spatial vertex domain.

Wavelet transform has been widely used in many fields such as signal processing, data compression, denoising, etc. Wavelet transform has the property that the wavelets are well localized in both frequency and space domain.

While Fourier transform has been successfully transfer to the graph domain ([Kipf & Welling 2017](https://arxiv.org/abs/1609.02907)), it is also possible to define wavelet transform on graphs. The spectral graph wavelet transform (SGWT) ([Hammond et al., 2011](https://www.sciencedirect.com/science/article/pii/S1063520310000552)) is one of the earliest works that brought the wavelet transform into the spectral graph theory.

The challenges of designing the wavelet transform on graphs are mainly twofold:

* **How to design the mother wavelet in graph domain, and how to define the scaling and translation operations for this mother wavelet on weighted graphs.**

  The mother wavelet is localized in some Euclidean space, however, graphs don't necessarily contains the coordinates information, then how can we design a localized mother wavelet in a non-Euclidean space?

  Another problem is, in the classical wavelet transform, both scaling and translation of the mother wavelet are defined in the space domain, however, for the graph data, it is elusive how to define these two operations directly.

  In the SGWT, the wavelet designing problem and the scaling problem are both solved by a wavelet transform analogue in the frequency domain. It is done by first showing that the scaling of a wavelet in the space domain can be completely transferred to the frequency domain, where a scaled wavelet acts as a bandpass filter (only frequencies in a frequency band are passed). Hence in the SGWT, instead of scaling in the space domain, we scale the mother wavelet in the frequency domain. A real-valued kernel function $g(\lambda)$ is designed to represent the corresponding bandpass filter of the mother wavelet in the frequency domain, hence we no longer need to worry about how to design the mother wavelet in a non-Euclidean space, we just need to make sure that the corresponding wavelets of $g(\lambda)$ in the graph domain are indeed localized. Then for some discrete scale values $s_1,s_2,\dots,s_J$, we could obtain multiple kernel functions $g(s_1\lambda), g(s_2\lambda),\dots,g(s_J\lambda)$ such that each of them will cover a sub-band of the entire graph spectrum.

  Then for the translation problem, a delta impulse function (only have value 1 at one point and 0 everywhere else) was introduced to show that translating a mother wavelet is the same as applying the mother wavelet to a delta impulse.

* **How to make the graph wavelet transform efficient for large graphs.**

  Based on the solution to the scaling problem, a Fourier transform of the graph data is needed, it is done by diagonalizing the graph Laplacian matrix $\mathscr{L}=U\Lambda U$ (also termed as the eigendecomposition) to obtain the frequencies of the graph and the corresponding eigenbasis. However, this eigendecomposition is of $O(N^3)$ time complexity. For any graph with more than a few thousands vertices, the SGWT will be too expansive to use.

  Hence, a polynomial approximation approach is used to calculate the transform without the expansive eigendecomposition. It is achieved by directly approximating the kernel function $g(s\lambda)$ by truncated Chebyshev polynomials. The maximum error of the truncated Chebyshev polynomials is only a slightly higher than that of the minimax polynomial (the unique polynomial which has the smallest maximum deviation from the true function), and in the region where $g(s\lambda)$ is smooth, truncated Chebyshev polynomials have significantly lower approximation error.

## Table of Contents
- [Classical Continuous Wavelet transform](#classical-continuous-wavelet-transform)
  * [CWT as a Bandpass Filter](#cwt-as-a-bandpass-filter)
  * [Delta Functions](#delta-functions)
- [Notation for Weighted Graphs](#notation-for-weighted-graphs)
- [Graph Fourier Transform](#graph-fourier-transform)
- [Spectral Graph Wavelet Transform (SGWT)](#spectral-graph-wavelet-transform--sgwt-)
  * [Scaling of Graph Wavelets](#scaling-of-graph-wavelets)
    + [Scaling in the Classical Continuous Wavelets Transform](#scaling-in-the-classical-continuous-wavelets-transform)
    + [Scaling Analogue on Weighted Graphs](#scaling-analogue-on-weighted-graphs)
  * [Translation of Graph Wavelets](#translation-of-graph-wavelets)
    + [Translation in the Classical Continuous Wavelets Transform](#translation-in-the-classical-continuous-wavelets-transform)
    + [Translation Analogue on Weighted Graphs](#translation-analogue-on-weighted-graphs)
- [Polynomial Approximation](#polynomial-approximation)
- [Inverse Graph Wavelet transform](#inverse-graph-wavelet-transform)
  * [Inverse of Continuous SGWT](#inverse-of-continuous-sgwt)
  * [Inverse of Discretized SGWT](#inverse-of-discretized-sgwt)
- [SGWT Kernel design](#sgwt-kernel-design)
- [Examples](#examples)
- [Appendix](#appendix)
  * [Frame Bound](#frame-bound)
  * [Spatial localization](#spatial-localization)
- [References](#references)

<small><i><a href='http://ecotrust-canada.github.io/markdown-toc/'>Table of contents generated with markdown-toc</a></i></small>

## Classical Continuous Wavelet transform
The classical continuous wavelet transform (CWT) is construct from the **mother wavelet function $\psi(x)$**, which is a wave-like function that will burst for a short time and then quickly die away. To generate a continuous family of the wavelets, we simply scale and translate the mother wavelet by a continuous scaling factor $s>0$ and a continuous translation factor $a$:

$$\psi_{s,a}(x)=\frac{1}{s}\psi(\frac{x-a}{s})$$

Note that the $\frac{1}{s}$ is the $L1$ normalization, a more common selection will be the $L2$ normalization ($\frac{1}{\sqrt{s}}$), the reason that $L1$ normalization is used in SGWT is because $L1$ normalization will make the wavelets at different scales all have the same amplitude in the frequency domain.

Then for given a signal $f(x)$, the **wavelet coefficients** are obtained by taken the inner product of each of these wavelet $\psi_{s,a}$ with the signal $f$:

$$W_f(s,a)=\langle \psi_{s,a}, f \rangle=\int_{-\infty}^{\infty}\frac{1}{s}\psi^*(\frac{x-a}{s})f(x)dx$$

We could think of the wavelet coefficient $W_f(s,a)$ at scale $s$ and translation $a$ being how much of $\psi_{s}$ exists in $f$ around $a$, with a continuous translation parameter $a$, we evaluate this amount of $\psi_{s}$ everywhere along $f$ with a fine step, and with a continuous scaling parameter $s$, we repeat this evaluation using a gradually changing wavelet $\psi_{s}$.

This process is invertible, i.e., $f(x)$ can be reconstructed using all the wavelet coefficients $W_f(s,a)$ at different scales and locations), provide that the mother wavelet $\psi(x)$ satisfies the **admissibility condition**:

$$\int_0^{\infty}\frac{|\hat{\psi}(\omega)|^2}{\omega}d\omega=C_\psi<\infty$$

Where $\hat{\psi}(\omega)$ being the Fourier transform of the mother wavelet $f(x)$. If this admissibility condition is satisfied, then we can invert the CWT:

$$f(x)=\frac{1}{C_\psi}\int_0^{\infty}\int_{-\infty}^{\infty}W_f(s,a)\psi_{s,a}(x)\frac{dads}{s}$$

So far, we are assuming that $f(x)$ is a signal on the real line, now we want to apply this classical CWT to graph data.

**The first challenge will be how to design a mother wavelet that is localized in the graph domain** as it is localized in space domain for the classical CWT, i.e., the graph wavelets should span around a center vertex in the graph with respect to some kind of distance measurement. Noted that this distance measurement should not depends on any kind of Euclidean distance since graph data normally don't have any Euclidean structure encoded. A naturally choice could be the shortest path between two vertices.

**The second challenge is how to define the scaling and translation operations for this mother wavelet on weighted graphs**. Assume that we are able to design a mother wavelet that is localized on a graph with respect to the shortest path distance measurement, then the scaling could be done by expanding/shrinking the radius of the mother wavelet, and the translation could be achieved by shifting the wavelet center to each vertex in the graph.

The SGWT did not approach to these two challenges directly, instead, they design the mother wavelet in the frequency domain, and shown that scaling could be also done in the frequency domain. As for the translation operation, a delta impulse function (only have value 1 at one point and 0 everywhere else) was introduced to show that, translating a mother wavelet is the same as applying the mother wavelet to a delta impulse, this is useful because we can easily design a delta function on graphs, which simply is a vector that only only have value 1 at certain vertex and 0 elsewhere.

Next we will first look at that, for the classical CWT, how can we transfer the scaling operation from the space domain to the frequency domain, and utilize delta impulse to achieved the translation of mother wavelet. Then in the following sections, after we got familiar with the notations of graph and the graph Fourier transform, we can discuss how to design the scaling and translation analogues on graphs.

### CWT as a Bandpass Filter
We already know that the CWT computes the inner product of a scaled and translated mother wavelet $\psi_{s,a}(x)=\frac{1}{s}\psi(\frac{x-a}{s})$ with a signal $f(x)$:

 $$CWT=\langle \psi_{s,a}, f \rangle=\int_{-\infty}^{\infty}\frac{1}{s}\psi^*(\frac{x-a}{s})f(x)dx$$

We can also interpret the CWT as a frequency-based filtering of the signal by rewriting the CWT as an inverse Fourier transform.

$$CWT=\frac{1}{2\pi}\int_{-\infty}^{\infty}\bar{\hat{\psi}}(s\omega)\hat{f}(\omega)e^{i\omega a}d\omega$$

From these two equations, we could see that stretching a wavelet in time will cause its support in the frequency domain to shrink, and its center frequency will move towards lower frequencies. The following figure demonstrates this effect for a hypothetical wavelet and scale factors of 1,2, and 4.

<div style="text-align: center">
<img src="https://www.mathworks.com/help/wavelet/gs/wavelet_freqdomain.png"/>
<p><i>Fig. 1. Fourier transform of the wavelets with scales 1,2 and 4. (Image source:<a href="https://www.mathworks.com/help/wavelet/gs/continuous-wavelet-transform-as-a-bandpass-filter.html">MathWorks</a>)</i></p>
</div>

In the above figure, $\omega$ corresponds to the frequency, $\hat{\psi}(\omega)$ is the Fourier transform of the wavelet function $\psi$, the amplitude of $\hat{\psi}(\omega)$ shows how much the wavelet contains certain frequency components, it is centered at $\omega_0$, which termed as the center frequency, it is the main frequency component of wavelet $\psi$. The support (none-zero part) of $\hat{\psi}(\omega)$ is the frequency band the wavelet contains (wavelet spectrum), when we do the wavelet transform, we multiple this $\hat{\psi}(\omega)$ with the frequencies of the signal ($\hat{f}$), only the frequency components within the support of the wavelet will remain, that is why CWT can act as a bandpass filter.

Also in this figure, from the top to the bottom are three wavelets with different scales (1, 2, 4 respectively), we see that the center frequency has changed along with the scale. Bigger the scale, lower the center frequency ($center frequency=\omega/s$). CWT coefficients at lower scales represent energy in the input signal at higher frequencies, while CWT coefficients at higher scales represent energy in the input signal at lower frequencies. Noticed that, the width of the bandpass filter also changed, it is inversely proportional to scale. The width of the CWT filters decreases with increasing scale. This follows from the uncertainty relationships between the time and frequency support of a signal: the broader the support of a signal in time, the narrower its support in frequency. The converse relationship also holds.

### Delta Functions
<div style="text-align: center">
<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/4/48/Dirac_distribution_PDF.svg/650px-Dirac_distribution_PDF.svg.png"/>
<p><i>Fig. 2. The delta function. (Image source:<a href="https://www.wikiwand.com/en/Dirac_delta_function#:~:text=As%20a%20distribution%2C%20the%20Dirac,of%20the%20Dirac%20delta%20function.">Wikipedia</a>)</i></p>
</div>

The **delta function $\delta(x)$**, also called the **Dirac delta function** or the **unit impulse symbol** in engineering and signal processing, is drawn as a singular arrow at one point to indicate an infinitely tall spike (an impulse). The y-axis value 1 represent its density ($\int_{-\infty}^{\infty}\delta(x)dx=1$), not height.

The delta function has an important property: **Convolution of a function $f(x)$ with a shifted delta function $\delta(x-a)$ yields a shifted version of that function of the same amount $f(x-a)$:**

$$f(x)\star\delta(x-a)=f(x-a)$$

Using this property, translating the mother wavelet by a factor $a$ is the same as convoluting the mother wavelet with a shifted delta function $\delta(x-a)$:

Let

$$\psi_{s}(x)=\frac{1}{s}\psi(\frac{x}{s})$$

then

$$\psi_{s}(x) \star \delta(x-a)=\psi_{s}(x-a)=\frac{1}{s}\psi(\frac{x-a}{s})=\psi_{s,a}(x)$$

Later, this result will be used to define the translation operation on graphs.

## Notation for Weighted Graphs
A weighted graph $G = `\{E, V , w\}`$ consists of:

* a set of vertices $V$, where $\|V\| = N < \infty$.
* a set of edges $E$.
* a weight function $w: E\rightarrow\mathbb{R}^+$ which assigns a positive weight to each edge.

The adjacency matrix $A$ for a weighted graph $G$ is the $N \times N$ matrix where

$$
A_{m,n} = \left
\{
\begin{array}{ll}
w(e), & \mbox{if $e\in E$ connects vertices $m$ and $n$} \\
0, & \mbox{otherwise} \\
\end{array}\right.
$$

The degree matrix $D$ is a $N \times N$ diagonal matrix Where

$$
D_{m,n} = \left
\{
\begin{array}{ll}
d(m)=\sum_nA_{m,n}, & \mbox{if $m=n$} \\
0, & \mbox{otherwise} \\
\end{array}\right.
$$

The Laplacian matrix, also called the graph Laplacian, is a matrix representation of a graph, many useful properties of a graph can be extract using the graph Laplacian. It is defined as

$$L=D-A$$

Noticed that:

* L is symmetric.
* L is positive-semidefinite (that is all its eigenvalues are non-negative).


## Graph Fourier Transform

## Spectral Graph Wavelet Transform (SGWT)
### Scaling of Graph Wavelets
#### Scaling in the Classical Continuous Wavelets Transform
#### Scaling Analogue on Weighted Graphs

### Translation of Graph Wavelets
#### Translation in the Classical Continuous Wavelets Transform
#### Translation Analogue on Weighted Graphs

## Polynomial Approximation

## Inverse Graph Wavelet transform
### Inverse of Continuous SGWT

### Inverse of Discretized SGWT

## SGWT Kernel design

## Examples

## Appendix
### Frame Bound

### Spatial localization

## References
[1] Hammond, David K., Pierre Vandergheynst, and Rémi Gribonval. "[Wavelets on graphs via spectral graph theory.](https://www.sciencedirect.com/science/article/pii/S1063520310000552)" Applied and Computational Harmonic Analysis 30.2 (2011): 129-150.

[2] Hammond, David K., Pierre Vandergheynst, and Rémi Gribonval. "[The spectral graph wavelet transform: Fundamental theory and fast computation.](https://hal.inria.fr/hal-01943589/document)" Vertex-Frequency Analysis of Graph Signals. Springer, Cham, 2019. 141-175.

[3] Kipf, Thomas N., and Max Welling. "[Semi-supervised classification with graph convolutional networks.](https://arxiv.org/abs/1609.02907)" arXiv preprint arXiv:1609.02907 (2016).
