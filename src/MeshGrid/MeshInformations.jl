abstract type AbstractMesh end

mutable struct Mesh <: AbstractMesh
  X::Matrix{Float64} # vector of nodes positions
  IEN::Matrix{Int64} # ID element -> ID nodes
  INE::Vector{Vector{Int64}} # ID node -> ID elements
  nsd::Int64 # number of spacial dimensions
  nnp::Int64 # number of all nodes
  nen::Int64 # number of element nodes
  nel::Int64 # number of all elements
  V_domain::Float64

  function Mesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
    V_domain = calculate_mesh_volume(X, IEN)

    X = reduce(hcat, X)
    IEN = reduce(hcat, IEN)
    INE = nodeToElementConnectivity(X, IEN)
    nsd = size(X, 1)
    nnp = size(X, 2)
    nen = size(IEN, 1)
    nel = size(IEN, 2)
   
    return new(X, IEN, INE, nsd, nnp, nen, nel, V_domain)
  end
end

function nodeToElementConnectivity(
  X::Matrix,
  IEN::Matrix,
)
  INE = [Vector{Int64}() for _ = 1:size(X, 2)]
  for el = 1:size(IEN, 2)
    for i = 1:size(IEN, 1)
      push!(INE[IEN[i, el]], el)
    end
  end
  return INE
end

mutable struct TriangularMesh <: AbstractMesh
  X::Matrix{Float64} # vector of nodes positions
  IEN::Matrix{Int64} # ID element -> ID nodes
  INE::Vector{Vector{Int64}} # ID node -> ID elements
  nsd::Int64 # number of spacial dimensions (3)
  nnp::Int64 # number of all nodes
  nen::Int64 # number of element nodes (3)
  nel::Int64 # number of all elements

  function TriangularMesh(X::Vector{Vector{Float64}}, IEN::Vector{Vector{Int64}})
    X = reduce(hcat, X)
    IEN = reduce(hcat, IEN)
    INE = nodeToElementConnectivity(X, IEN)
    nsd = size(X, 1)
    nnp = size(X, 2)
    nen = size(IEN, 1)
    nel = size(IEN, 2)

    return new(X, IEN, INE, nsd, nnp, nen, nel)
  end
end
