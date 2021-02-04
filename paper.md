title: 'kallisto: A Python command-line interface for modelling in computational chemistry and atomic features in machine learning applications'
tags:
  - Python
  - Poetry
  - Click
  - Chemistry
  - Machine Learning
  - Atomic Features
  - Computational Chemistry
  - Computational Modelling
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

Machine learning (ML) has recently become even more popular within the pharmaceutical industry especially for tasks as, e.g., building predictive models, performing virtual screening, or predicting compound activities.[`@ghasemi2018`]
Traditional ML models often exploited the quantitative structure-activity relationship (QSAR) that has been popularized by medicinal chemists and statisticians to relate bioactivities to specific functional group manipulations.
The QSAR approach decreases the dimensionalty of the underlying problem and projects the molecular structure into a vector space spanned by the physicochemical features.
Early approaches heavily relied on linear regression to correlate physicochemical features and bioactivities for each compound.
Modern approaches, however, combine such features with non-linear machine learning algorithms and driving new scales of chemical research.

There exist great chemoinformatic packages like RDKit[`@landrum2006`] that enable the calculation of physicochemical features based on chemoinformatical information like the molecular graph.
However, often we want to go beyond a structure-only approach to incorporate electronic structure effects that are described by a quantum mechanical treatment.
The calculation of such kind of physicochemical features often applies well-established quantum chemistry methods like Kohn-Sham density functional theory (DFT) that is currently the workhorse of computational chemistry.[`@parr1980`, `@kohn1999`]
However, generating the feature space by DFT is computationally demanding and can become the computational bottleneck especially when aiming for high-throughput experiments with several hundrets to thousands of molecules.
Since there exists a critical need for an efficient featurizer, we developed the ``kallisto`` command-line interface that solves this problem by enabling a fast way to generate various types of physicochemical atomic features for atoms and molecules for most of the periodic system.
Molecular geometries need to have an [``xmol``](https://en.wikipedia.org/wiki/XYZ_file_format) or a [Turbomole](https://www.turbomole.org/wp-content/uploads/2019/11/Turbomole_Manual_7-4-1.pdf) like format to be processed by ``kallisto``.
In addition we implemented several computational modelling helper to simplify the creation of high-throughput procedures.
Some of those modelling helper depend on the open source [xtb](https://github.com/grimme-lab/xtb) tight binding scheme developed by Stefan Grimme and co-worker.[`@bannwarth2020`]
The [online documentation](https://ehjc.gitbook.io/kallisto/) covers all high-level features of this software mostly in terms of copy-paste recipes.
Furthermore, we cover parts of the theory that are used behind the scenes and compare to experimental data as well as to modern deep learning models.


# Acknowledgements

EC acknowledges contributions from Philipp Pracht (`@ppracht`).

# References
