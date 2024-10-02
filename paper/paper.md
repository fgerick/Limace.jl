---
title: 'Limace.jl: A Julia package to compute hydromagnetic modes in spherical domains'
tags:
  - Julia
  - geophysics
  - fluid dynamics
  - geomagnetism
  - hydromagnetic modes
authors:
  - name: Felix Gerick
    orcid: 0000-0001-9924-0562
    equal-contrib: true
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
affiliations:
 - name: National Centre for Space Studies, France
   index: 1
   ror: 04h1h0y33
 - name: Royal Observatory of Belgium, Belgium
   index: 2
   ror: 00hjks330
date: 2 October 2024
bibliography: paper.bib
---

# Summary

[@gerickinterannual2024]

# Statement of need


# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Figures

Figures can be included like this:
![Torsional Alfven mode.\label{fig:tm}](torsionalmodes.png){ width=50%}
and referenced from text using \autoref{fig:example}.

# Acknowledgements

I have received funding from the European Research Council (ERC) GRACEFUL Synergy Grant No. 855677. 
This project has been funded by ESA in the framework of EO Science for Society, through contract 4000127193/19/NL/IA (SWARM + 4D Deep Earth: Core). 
I thank Phil Livermore for the help in the theoretical development of the model.

# References
