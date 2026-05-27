const Dimensionable = Union{AbstractMove, AbstractLattice}

dim(dd::Dimensionable) = dd.d
dim(dd::AbstractPoints) = length(dd.points[1])
dim(em::EstimationMethod) = dim(em.pl)

#dim(h::HamiltonianFamily) = dim(h.hs[1])
dim(h::AbstractHamiltonian) = h isa HamiltonianFamily ? dim(h.hs[1]) : dim(h.pointset)
