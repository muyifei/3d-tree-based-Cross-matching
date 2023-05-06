# 3d-tree-based-Cross-matching
an optimized cross-matching algorithm using 3d-tree

# abstract
In this paper, we propose a 3d-tree-based cross-matching algorithm that converts the angular distance formula into an equivalent 3d Euclidean distance and uses 3d-tree method to reduce the overall computational complexity and to avoid boundary issues. Furthermore, we demonstrate the superiority of the 3d-tree approach over the 2d-tree method and implement it using a multi-threading technique during both the construction and querying phases. We have experimentally evaluated the proposed 3d-tree cross-matching algorithm using publicly available catalog data. The results show that our algorithm applied on two 32-core CPU achieves equivalent performance than previous experiments conducted on a six-node CPU-GPU cluster.

# tutorial
Firstly, we need two catalog data files to cross-match. Each file involves right ascension and declination information within binary format of c++ double(8 bytes), and each record holds 16 bytes. The input data should fulfill this constraint: ra in [0-360) and dec in [-90-90]
