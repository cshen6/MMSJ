Wikipedia Documents Matlab Data
=====================================================

The Matlab data is made available for its use together with other Matlab code in the following paper: 

Cencheng Shen and Carey E. Priebe, "Manifold Matching using Shortest-Path Distance and Joint Neighborhood Selection", submitted, on arxiv, http://arxiv.org/abs/1412.4098, 2015.

Please feel free to use the data for non-commercial purposes, as long as proper acknowledgement and citation are made. 

For any feedback, please contact cshen6@jhu.edu, youngser@jhu.edu, or cep@jhu.edu.

Contributors in Data Collection & Processing
(in alphabetical order of last name): 
Dave Marchette, Youngser Park, Carey Priebe, Cencheng Shen, Ming Sun

Affiliations:
Department of Applied Mathematics and Statistics & Center of Imaging Science,
Johns Hopkins University,  
Baltimore, Maryland, United States.



Information
-------------------------

The Wikipedia.mat data collects text and network features for n=1382 documents, which are within 2-neighborhood of the English article "algebraic geometry" back in 2011. The documents are labelled into 5 classes.

It contains the following matrix: TE, TF, GE, GF, GEAdj, GFAdj, Label.

TE is the text feature distance matrix of the documents in English. After stopwords removal and words stemming, the latent semantic indexing feature is calculated for each document, followed by cosine distance to yield the distance matrix.

TF is the text feature distance matrix of those documents in French.

GEAdj contains the adjacency information of those documents in Wikipedia English.

GFAdj contains the adjacency information of those documents in Wikipedia French.

GE processes GEAdj into shortest-path distances, with any path longer than 4 imputed to 6.

GF processes GFAdj into shortest-path distances, with any path longer than 4 imputed to 6.

Label contains the label information for the documents, which has the following classes:
0 -- Category
1 -- People
2 -- Locations
3 -- Dates
4 -- Math
The classification is manually done to the best of our knowledge.



Relevant Papers
-------------------------

We provide a complete list of our papers that have used this data. And we will appreciate that if you are willing to cite any of them.

C. Shen and C. E. Priebe, ※Manifold Matching using Shortest-Path Distance and Joint Neighborhood Selection,§ submitted, 2015.

C. Shen, M. Sun, M. Tang, and C. E. Priebe, ※Generalized canonical correlation analysis for classification,§ Journal of Multivariate Analysis, vol. 130, pp. 310每322, 2014.

M. Sun and C. E. Priebe, ※Efficiency investigation of manifold matching for text document classification,§ Pattern Recognition Letters, vol. 34, no. 11, pp. 1263每1269, 2013.

M. Sun, C. E. Priebe, and M. Tang, ※Generalized canonical correlation analysis for disparate data fusion,§ Pattern Recognition Letters, vol. 34, no. 2, pp. 194每200, 2013.

C. E. Priebe, D. J. Marchette, Z. Ma, and S. Adali, ※Manifold matching: Joint optimization of fidelity and commensurability,§ Brazilian Journal of Probability and Statistics, vol. 27, no. 3, pp. 377每400, 2013.

