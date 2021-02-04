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

Machine learning (ML) has recently become very popular within pharmaceutical industry.
Especially tasks as, e.g., building predictive models, performing virtual screening, or predicting compound activities are great usecases for potential ML applications.[`@roy2015`, `@sprous2010`]
Traditional ML models often rely on the quantitative structure-activity relationship (QSAR) that has been popularized by medicinal chemists and statisticians to relate bioactivities to specific functional group manipulations.
This QSAR approach decreases the dimensionality of the underlying problem and projects the molecular structure into a vector space spanned by the physicochemical features.
While early approaches relied more on linear regression, modern approaches combine features with non-linear machine learning algorithms.

Great chemoinformatic packages like RDKit[`@landrum2006`] enable the calculation of physicochemical features based on structural information like the molecular graph.
However, often we want to go beyond a structure-only approach thus incorporating electronic structure effects as obtained, e.g., by a quantum mechanical (QM) treatment.
The calculation of QM-derived features relies often on well-established quantum chemistry methods like Kohn-Sham density functional theory (DFT) that is currently the workhorse of computational chemistry.[`@parr1980`, `@kohn1999`]
However, generating the feature space by DFT is computationally demanding and can become the computational bottleneck especially when aiming for high-throughput experiments with several hundred to thousands of molecules.

Since there exists a critical need for an efficient yet accurate atomic featurizer, we developed the ``kallisto`` command-line interface that is able to calculate QM-derived atomic features for atoms and molecules efficiently (whole periodic table up to Radon).
Molecular geometries need to have an [``xmol``](https://en.wikipedia.org/wiki/XYZ_file_format) or a [``Turbomole``](https://www.turbomole.org/wp-content/uploads/2019/11/Turbomole_Manual_7-4-1.pdf) like format to be processed by ``kallisto``.
Besides, we implemented several computational modelling helpers to simplify the creation of high-throughput procedures.
Some of those modelling helpers depend on the open-source [xtb](https://github.com/grimme-lab/xtb) tight-binding scheme that has been developed by Stefan Grimme and co-worker.[`@bannwarth2020`]
The [online documentation](https://ehjc.gitbook.io/kallisto/) covers all high-level functionalizations of this software mostly in terms of copy-paste recipes.
Furthermore, we cover bits of the underlying theory and compare to experimental data as well as to other modern deep learning models.


# Acknowledgements

EC acknowledges contributions from Philipp Pracht (`@pprcht`).

# References
