# Structural-Polarization-Twitter

## Final Research Project - Complex and Social Networks (CSN)

**Authors:** Sara Montese, Marius Behret 

## Table of Contents
1. [Domain](#domain-background)
2. [Aims of the Research](#aims-of-the-research)
3. [Hypotheses to be Tested](#hypotheses-to-be-tested)
4. [Theoretical Framework](#theoretical-framework)
5. [Methods](#methods)
6. [Report](#report)
7. [Code](#code)


![rt_israel edges](https://github.com/SaraMo14/Structural-Polarization-Twitter/assets/74814020/45dd1779-56d5-4f3d-9ddc-fb26d6e30990)

## Domain
Structural polarization refers to the phenomenon where a network (such as a social or political network) becomes divided into distinct groups or “poles” based on certain opinions. 
In the context of the paper “Separating Polarization from Noise: Comparison and Normalization of Structural Polarization Measures” [1], structural polarization refers to the quantification of polarization within a network structure.

## Aims of the Research
 - Analyze structural polarization of datasets from Social Media
 - Compare the results of standard and normalized polarization measures
 - Study the correlation of the polarization measures under analysis
 - Propose one or more polarization scores

## Hypotheses to be Tested
- Randomized networks appear to be polarized even if they are randomly generated using the eight polarization measures from [1]
- The choice of the polarization scores has a lower impact than their normalization
- The process of normalizing polarization scores effectively reduces significant interference caused by the local properties of the network

## Theoretical Framework
Network science, social network analysis, polarization measurement, community detection, and clustering.

## Methods
Generation of the randomized non-polarized networks; calculation of the several polarization scores like Random Walk Controversy (RWC), Adaptive Random Walk Controversy (ARWC), Betweenness Centrality Controversy (BCC), Boundary Polarization (BP), Dipole polarization (DP), E-I Index (EI), Adaptive E-I Index (AEI) and Modularity (Q) on polarized and non-polarized networks. Comparison between standard and normalized polarization metrics.

## Report
For more detailed information, please refer to our research report [here](./report/CSN_Final_Project.pdf).

## Code
Find the code implementation in the 'code' directory.

## References
1. A. Salloum, T. H. Y. Chen, and M. Kivelä, “Separating polarization from noise: comparison and normalization of structural polarization measures,” Proceedings of the ACM on human-computer interaction, vol. 6, no. CSCW1, pp. 1–33, 2022

## Acknowledgments
This research is conducted at the Facultat d'Informàtica de Barcelona (FIB), Universitat Politècnica de Catalunya (UPC) - BarcelonaTech.


