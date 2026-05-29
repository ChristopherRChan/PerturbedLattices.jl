const Dimensionable = Union{AbstractMove, AbstractLattice}

dim(dd::Dimensionable) = dd.d
dim(dd::AbstractPoints) = length(dd.points[1])
dim(dd::AbstractPerturbedLatticeModel) = length(dd.pointset.points[1]) 
dim(em::EstimationMethod) = dim(em.pl)