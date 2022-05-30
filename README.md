# dS_Bell
Study of Bell inequalities in real space for de Sitter spacetime

This is supplementary code to the paper "Real-space Bell inequalities in de Sitter", a research project developed in collaboration with Vincent Vennin from APC Paris. The paper is under review by JCAP, but you can already check it on the arXiv https://arxiv.org/abs/2203.03505. The abstract goes like this

"Bell-inequality violations reveal the presence of quantum correlations between two particles that have interacted and then separated. Their generalisation to quantum fields is necessary to study a number of field-theoretic setups, such as cosmological density fluctuations. In this work, we show how Bell operators can be constructed for quantum fields in real space, and for Gaussian states we compute their expectation value in terms of the field power spectra. We then apply our formalism to a scalar field in de-Sitter space-time. We find that, in spite of the tremendous production of entangled particles with opposite wave momenta on large scales, Bell inequalities are not violated in real space. The reason is that, when considering measurements of a field at two distinct locations in real space, one implicitly traces over the configuration of the field at every other location, leading to a mixed bipartite system. This "effective decoherence" effect is responsible for the erasure of quantum features, and casts some doubts on our ability to reveal the quantum origin of cosmological structures. We finally discuss these results in the light of quantum discord."

This repository contains the code used to produce the plots displayed in the paper. As such, it can

1) Compute the covariance matrix of the quantum state of a scalar field in Minkowski and de Sitter spacetimes.
2) Use this covariance matrix to compute the expected value of Bell operators based on pseudo-spin operators.
3) Test Bell inequalities in real-space (i.e., position space) for cosmological perturbations.
