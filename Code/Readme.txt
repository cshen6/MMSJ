Manifold Matching Matlab Code
=====================================================

Author: Cencheng Shen

Affiliations: 
Department of Applied Mathematics and Statistics & Center of Imaging Science,
Johns Hopkins University,  
Baltimore, Maryland, United States.

Date: 1/13/2015

Version: 0.1



Copyrights
-------------------------
(C) Cencheng Shen, 2015

The code is for the following paper: 

Cencheng Shen and Carey E. Priebe, "Manifold Matching using Shortest-Path Distance and Joint Neighborhood Selection", submitted, on arxiv, http://arxiv.org/abs/1412.4098, 2015.

Feel free to modify, extend or distribute this code for non-commercial purposes, as long as proper acknowledgement and citation are made. 

And please contact cshen6@jhu.edu or cep@jhu.edu for any feedback and bug reports.



Installation
-------------------------
Add the entire folder into Matlab running path (add with subfolders), and the functions can be used.



Functionality
-------------------------
The main purpose of the code is to do manifold matching based on our paper, and the main functions are ManifoldMatching.m and ManifoldMatchingEuc.m.  

Within them, neighborhoodIsomap.m and neighborhoodLLE.m do Isomap and LLE for distance matrices. GeneralCCA.m does generalized canonical correlation analysis based on iterative minimization of sum of pairwise correlation. GeneralProcrustes.m finds the Procrustes matchings using greedy approximation. SMDS.m does raw-sress multidimensional scaling, with OOSMDS.m doing out-of-sample embedding for MDS. And we also implement testing power for matching (plotPower.m), with other auxliary functions (plotVelocity.m, GetRealData.m, GetRealDataEuc.m) useful for visualization and data formatting.

We also include an amazing dimension reduction toolbox from L. Maaten, and we delete the outlier detection step in the toolbox for our matching purpose. 

The running examples on http://www.cis.jhu.edu/~cshen/ provides some illustrations on how to use the code.



Acknowledgements
-------------------------
We have used the Isomap code from J. Tenenbaum, and the LLE code from S. Roweis and L. Saul, with heavy modifications for our purpose (mainly adding joint neighborhood selection and out-of-sample extensions). 

And LTSA/Laplacian eigenmaps/HLLE functionality is based on the dimension reduction toolbox of L. Maaten.

Links to their websites are also provided on our download webpage.


