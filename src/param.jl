# params

const Parametrisable = Union{AbstractPerturbedLatticeModel, AbstractMove, AbstractHamiltonian, EstimationMethod}

function params end
function nbparam end
function params! end

# alias
θ = params
nθ = nbparam
θ! = params!

