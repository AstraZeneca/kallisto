---
title: 'kallisto: A command-line interface to simplify computational modelling and the generation of atomic features.'
tags:
  - Python
  - Poetry
  - Atomic Features
  - Machine Learning
  - Computational Chemistry
  - Pharmaceutical Science

authors:
  - name: Eike Caldeweyher^[Corresponding author]
    orcid: 0000-0002-3985-595X
    affiliation: "1" 

affiliations:
 - name: Data Science and Modelling, Pharmaceutical Science, R&D, AstraZeneca, Gothenburg, Sweden
   index: 1

date: 4 February 2021

bibliography: paper.bib

---

# Statement of Need

Machine learning (ML) has recently become very popular within pharmaceutical industry [@roy2015] [@sprous2010].
Tasks as, e.g., building predictive models, performing virtual screening, or predicting compound activities are potential usecases for such ML applications.
Traditionally, ML models often rely on the quantitative structure-activity relationship (QSAR) that has been popularized by medicinal chemists and statisticians to relate bioactivities to specific functional group manipulations [@dudek2006] [@verma2010].
This QSAR approach decreases the dimensionality of the underlying problem and projects the molecular structure into a space spanned by the physicochemical features.
While early approaches relied more on linear regression, modern approaches combine such features with non-linear ML algorithms.

Chemoinformatic packages like RDKit [@landrum2006] enable the fast calculation of atomic/molecular features based on structural information like the molecular graph.
However, frequently we want to go beyond a structure-only approach thus incorporating electronic structure effects as obtained, e.g., by a quantum mechanical (QM) treatment.
The calculation of QM-based features relies often on well-established quantum chemistry methods like Kohn-Sham density functional theory (DFT) that is currently the workhorse of computational chemistry [@parr1980] [@kohn1999].
However, generating the feature space by DFT is computationally demanding and can become the computational bottleneck especially when aiming for high-throughput experiments with several hundred to thousands of molecules.

Since there exists a critical need for an efficient yet accurate featurizer, we developed the ``kallisto`` command-line interface that is able to calculate QM-based atomic features for atoms and molecules efficiently (whole periodic table up to Radon).
The features are either interpolating high-level time-dependent DFT [@dreuw2005] data (e.g., static/dynamic polarizabilitities with TD-PBE0 [@adamo1999] / d-aug-def2-QZVP [@weigend2005] references) or are parametrized [@caldeweyher2019] to reproduce QM-data (e.g., atomic partial charges with PBE0 [@adamo1999] / def2-TZVP [@weigend2005] Hirshfeld [@hirshfeld1977] references).
Molecular geometries need to have an [``xmol``](https://en.wikipedia.org/wiki/XYZ_file_format) or a [``Turbomole``](https://www.turbomole.org/wp-content/uploads/2019/11/Turbomole_Manual_7-4-1.pdf) like format to be processed by ``kallisto``.
Besides, we implemented several computational modelling helpers to simplify the development of high-throughput procedures.
Some of those modelling helpers depend on the open-source [xtb](https://github.com/grimme-lab/xtb) tight-binding scheme that has been developed by Stefan Grimme and co-wo
rker [@bannwarth2020].
The [online documentation](https://ehjc.gitbook.io/kallisto/) covers all high-level functionalizations of this software mostly in terms of copy-paste recipes.
Furthermore, we cover bits of the underlying theory and compare to experimental data as well as to other modern deep learning models.

# Acknowledgements

EC acknowledges contributions from Philipp Pracht (`@pprcht`).

# References
