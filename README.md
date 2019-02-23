# Spin Transport Research in Graphene

This is a fortran program running on Linux clusters (typically with 24 nodes) to solve spin transport problem with scattering and spin-orbit coupling in graphene.

This is a hybrid parallel computation program. It applies MPI to perform distributed parallel computation on cluster nodes and OPENMP to perform memory-shared parallel computation on each node.
