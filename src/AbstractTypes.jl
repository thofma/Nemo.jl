###############################################################################
#
#   Abstract types
#
###############################################################################

abstract Pari

abstract Flint

abstract Antic

abstract Generic

abstract Set{T}

abstract Ring{T} <: Set{T}

abstract Field{T} <: Ring{T}

abstract SetElem

abstract RingElem <: SetElem

abstract FieldElem <: RingElem

abstract PolyElem{T} <: RingElem

abstract ResidueElem{T} <: RingElem

abstract FractionElem{T} <: FieldElem

abstract SeriesElem{T} <: RingElem

# not always mathematical ring elements
abstract MatElem <: RingElem

abstract FiniteFieldElem <: FieldElem

abstract NumberFieldElem <: FieldElem

abstract MaximalOrderElem <: RingElem
