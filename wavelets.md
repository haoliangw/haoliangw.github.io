# Spectral Graph Wavelet Transform

> The spectral graph wavelet transform (SGWT) ([Hammond et al., 2011](https://www.sciencedirect.com/science/article/pii/S1063520310000552)) defines the wavelet transform for weighted graphs. The key point of the SGWT is to design a real-valued kernel function $g(\lambda)$ that acts as a bandpass filter in the Fourier domain. With the kernel functions $g(s\lambda)$ at discrete scales $s$, we are able to create wavelets that are not only localized in the graph spectrum, but also localized in the spatial vertex domain.

Wavelet transform has been widely used in many fields such as signal processing, data compression, denoising, etc. Wavelet transform has the property that the wavelets are well localized in both frequency and space domain.

While Fourier transform has been successfully transfer to the graph domain ([Kipf & Welling 2017](https://arxiv.org/abs/1609.02907)), it is also possible to define wavelet transform on graphs. The spectral graph wavelet transform (SGWT) ([Hammond et al., 2011](https://www.sciencedirect.com/science/article/pii/S1063520310000552)) is one of the earliest works that brought the wavelet transform into the spectral graph theory.

There are two major challenges of designing the wavelet transform on graphs:

* **How to define the scaling and translation operations of the mother wavelet on graphs.**

  In the classical wavelet transform, both scaling and translation of the mother wavelet are defined in the spatial domain, however for the graph data, it is elusive how to define these two operations directly.

  In the SGWT, this is solved by first showing that the scaling of a wavelet in the space domain can be completely transferred to the frequency domain, where a scaled wavelet acts as a bandpass filter (only frequencies in a frequency band are passed). Hence in the SGWT, instead of scaling in the space domain, we scale the mother wavelet in the frequency domain. A real-valued kernel function $g(\lambda)$ is designed to represent the corresponding bandpass filter of the mother wavelet in the frequency domain, then for some discrete scale values $s_1,s_2,\dots,s_J$, we could obtain multiple kernel functions $g(s_1\lambda), g(s_2\lambda),\dots,g(s_J\lambda)$ such that each of them will cover a sub-band of the entire graph spectrum.

  Then for the translation operation, a delta impulse function (only have value 1 at one point and 0 everywhere else) is used to show that translation of the mother wavelet is the same as applying the mother wavelet to a delta impulse.

* **How to make the graph wavelet transform efficient for large graphs.**

  The SGWT is based on the eigendecomposition of a graph Laplacian matrix $\mathscr{L}=U\Lambda U$. However, this eigendecomposition is of $O(N^3)$ time complexity. For any graph with more than a few thousands vertices, the SGWT will be too expansive to use.

  Hence, a polynomial approximation approach is used to calculate the transform without the expansive eigendecomposition. It is achieved by directly approximating the kernel function $g(s\lambda)$ by truncated Chebyshev polynomials. The maximum error of the truncated Chebyshev polynomials is only a slightly higher than that of the minimax polynomial (the unique polynomial which has the smallest maximum deviation from the true function), and in the region where $g(s\lambda)$ is smooth, truncated Chebyshev polynomials have significantly lower approximation error.

## Table of Contents
- [Classical Continuous Wavelet transform](#classical-continuous-wavelet-transform)
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

## Notation for Weighted Graphs

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
