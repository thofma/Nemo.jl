using Nemo

include("Polygons.jl")

################################################################################
#
#   Types
#
################################################################################

# Type for the parent ring
type PuiseuxSeriesRing{T} <: Nemo.Ring
  base_ring::Nemo.Ring
  S::Symbol
end

function PuiseuxSeriesRing(R, s)
  S = PuiseuxSeriesRing{elem_type(R)}(R, Symbol(s))
  return S, gen(S)
end

# Type for the Puiseux series 
abstract AbsPuiseuxSeries <: Nemo.RingElem

type PuiseuxSeries{T} <: AbsPuiseuxSeries
  coeff::Array{T, 1}
  exp::Array{Rational{Int}, 1}
  parent::Nemo.Ring
  prec::Rational{Int}

  function PuiseuxSeries{T}(coeff::Array{T, 1}, exp::Array{Rational{Int}, 1})
    z = new{T}()
    z.coeff = coeff
    z.exp = exp
    return z
  end
end

################################################################################
#
#   Basic access
#
################################################################################

Nemo.base_ring(R::PuiseuxSeriesRing) = R.base_ring

Nemo.parent(a::PuiseuxSeries) = a.parent

################################################################################
#
#   Functions for turning it into a proper ring for Nemo
#
################################################################################

Nemo.parent_type{T}(::Type{PuiseuxSeries{T}}) = PuiseuxSeriesRing{T}

Nemo.needs_parentheses(::PuiseuxSeries) = true

Nemo.show_minus_one{T}(::Type{PuiseuxSeries{T}}) = true

Nemo.is_negative(::PuiseuxSeries) = false

Nemo.elem_type{T}(::PuiseuxSeriesRing{T}) = PuiseuxSeries{T}

Nemo.zero{T}(R::PuiseuxSeriesRing{T}) = R(T[], Rational{Int}[])

Nemo.one(R::PuiseuxSeriesRing) = R(1)

(R::PuiseuxSeriesRing)() = zero(R)

Nemo.iszero(a::PuiseuxSeries) = (a == parent(a)(0))

Nemo.isone(a::PuiseuxSeries) = (a == parent(a)(1))

################################################################################
#
#   Constructors
#
################################################################################

function (a::PuiseuxSeriesRing{T}){T}(coeff::Array{T, 1}, exp::Array{Rational{Int}, 1})
  z = PuiseuxSeries{T}(coeff, exp)
  z.parent = a
  return z
end

function (a::PuiseuxSeriesRing{T}){T}(coeff::Array{T, 1}, exp::Array{Int, 1})
  return a(coeff, convert(Array{Rational{Int}, 1}, exp))
end

(a::PuiseuxSeriesRing){S <: Integer}(coeff::Array{S, 1}, exp) = a(map(a.base_ring, coeff), exp)

(a::PuiseuxSeriesRing)(b::Integer) = b != 0 ? a([ a.base_ring(b)] , [ 0//1 ]) : zero(a)

(a::PuiseuxSeriesRing{T}){T}(b::T) = a([b], [0//1])

Nemo.gen(R::PuiseuxSeriesRing) = R([R.base_ring(1)], [1//1])

################################################################################
#
#   Addition
#
################################################################################

function Base.:+{T}(a::PuiseuxSeries{T}, b::PuiseuxSeries{T})
  new_coeff = Array{T}(0)
  new_exp = Array{Rational{Int}}(0)

  ia, ib = 1, 1;

  p = min(a.prec, b.prec)

  # Append finite terms of a and b in proper order (=merge sort of exponent arrays) 
  # until min(prec(a), prec(b)).
  while ia <= length(a.coeff) && ib <= length(b.coeff)
    if a.exp[ia] < b.exp[ib]
      push!(new_exp, a.exp[ia])
      push!(new_coeff, a.coeff[ia])
      ia = ia + 1
    elseif a.exp[ia] > b.exp[ib]
      push!(new_exp, b.exp[ib])
      push!(new_coeff, b.coeff[ib])
      ib = ib + 1
    else
      c = a.coeff[ia] + b.coeff[ib]
      if !iszero(c)
        push!(new_exp, a.exp[ia])
        push!(new_coeff, c)
      end
      ia = ia + 1
      ib = ib + 1
    end
  end

  # If prec(b) < prec(a), add finite terms of a below prec(b)
  for i in ia:length(a.coeff)
    if !inprecision(b, a.exp[i])
      break
    end
    push!(new_exp, a.exp[i])
    push!(new_coeff, a.coeff[i])
  end

  # If prec(a) < prec(b), add finite terms of b below prec(a)
  for i in ib:length(b.coeff)
    if !inprecision(a, b.exp[i])
      break
    end
    push!(new_exp, b.exp[i])
    push!(new_coeff, b.coeff[i])
  end

  # Determine precision of a + b as the minimum (if it exists)
  if finiteprecision(a)
    if finiteprecision(b)
      push!(new_exp, min(a.exp[end], b.exp[end]))
    else
      push!(new_exp, a.exp[end])
    end
  else
    if finiteprecision(b)
      push!(new_exp, b.exp[end])
    end
  end

  return parent(a)(new_coeff, new_exp)
end

################################################################################
#
#   Multiplication
#
################################################################################

function Base.:*{T}(a::PuiseuxSeries{T}, b::PuiseuxSeries{T})
  if length(a.coeff) == 0 || length(b.coeff) == 0
    return zero(parent(a))
  end

  # Determine precision p of a*b,
  # p > highest exponent(a) + highest exponent(b) means precision of a*b is infinite 
  if finiteprecision(a)
    if finiteprecision(b)
      p = min(valuation(a) + precision(b), valuation(b) + precision(a))
    else
      p = valuation(b) + precision(a)
    end
  else
    if finiteprecision(b)
      p = valuation(a) + precision(b)
    else
      p = a.exp[end] + b.exp[end] + 1
    end
  end

  # Compute the product of all terms (assuming sum of exponent less than p) 
  # and sort the result by their exponent
  new_exp = Array{Rational{Int}}(0)
  new_coeff = Array{T}(0)
  for i_a in 1:length(a.exp)
    for i_b in 1:length(b.exp)
      e = a.exp[i_a] + b.exp[i_b]
      if e < p
        push!(new_exp, e)
        c = a.coeff[i_a] * b.coeff[i_b]
        push!(new_coeff, c)
      end
    end
  end
  v = sortperm(new_exp)
  
  # construct a*b, merging all coefficients with the same exponent 
  new_new_exp = Array{Rational{Int}}(0)
  new_new_coeff = Array{T}(0)

  push!(new_new_exp, new_exp[v[1]])
  push!(new_new_coeff, new_coeff[v[1]])

  for i in 2:length(v)
    if new_exp[v[i]] == new_new_exp[end]
      new_new_coeff[end] += new_coeff[v[i]]
    else
      if iszero(new_new_coeff[end])
        new_new_exp[end] = new_exp[v[i]]
        new_new_coeff[end] = new_coeff[v[i]]
      else
        push!(new_new_exp, new_exp[v[i]])
        push!(new_new_coeff, new_coeff[v[i]])
      end
    end
  end

  if iszero(new_new_coeff[end])
    pop!(new_new_exp)
    pop!(new_new_coeff)
  end

  # if either a or b was finite precision, set result to be of finite precision p
  if finiteprecision(a) || finiteprecision(b)
    push!(new_new_exp, p)
  end

  return parent(a)(new_new_coeff, new_new_exp)
end

################################################################################
#
#   Unsafe operations
#
################################################################################

function Nemo.mul!(a::PuiseuxSeries, b::PuiseuxSeries, c::PuiseuxSeries)
  z = b * c
  a.exp = z.exp
  a.coeff = z.coeff
  return nothing
end

function Nemo.add!(a::PuiseuxSeries, b::PuiseuxSeries, c::PuiseuxSeries)
  z = b + c
  a.exp = z.exp
  a.coeff = z.coeff
  return nothing
end

function Nemo.addeq!(a::PuiseuxSeries, b::PuiseuxSeries)
  z = a + b
  a.exp = z.exp
  a.coeff = z.coeff
  return nothing
end

################################################################################
#
#   Exponentiation
#
################################################################################

function Base.:(^)(a::PuiseuxSeries, b::Rational{Int})
  a != gen(parent(a)) && error("Not supported")

  return parent(a)([one(parent(a).base_ring)], [b])
end

################################################################################
#
#   Ad hoc operations
#
################################################################################

function Base.:(*){T <: RingElem}(a::PuiseuxSeries{T}, b::T)
#  if iszero(b)
#    return zero(parent(a))
#  end
  return parent(a)([ b * c for c in a.coeff], a.exp)
end

Base.:(*){T <: RingElem}(a::T, b::PuiseuxSeries{T}) = b * a

Base.:(*)(a::Integer, b::PuiseuxSeries) = parent(b).base_ring(a)*b

Base.:(*)(a::PuiseuxSeries, b::Integer) = b * a

################################################################################
#
#   Comparison
#
################################################################################

function Base.:(==)(a::PuiseuxSeries, b::PuiseuxSeries)
  return a.exp == b.exp && a.coeff == b.coeff
end

Base.:(==)(a::PuiseuxSeries, b::Integer) = (a == parent(a)(b))

Base.:(==)(a::Integer, b::PuiseuxSeries) = (b == a)

################################################################################
#
#  Random stuff
#
################################################################################

valuation(a) = a.exp[1]

finiteprecision(a) = length(a.coeff) != length(a.exp)

precision(a) = a.exp[end]

function inprecision(a, e) 
  if !finiteprecision(a)
    return true
  else
    return e < a.exp[end]
  end
end

function initial(a)
  iszero(a) && error("a must be nonzero: $a")
  return a.coeff[1]
end

################################################################################
#
#   String I/O
#
################################################################################

function Base.show(io::IO, R::PuiseuxSeriesRing)
  print(io, "Puiseux series ring in $(string(R.S)) over $(R.base_ring)")
end

# todo, don't print 1 in 1 * t
function Base.show(io::IO, f::PuiseuxSeries)
  t = parent(f).S

  if length(f.coeff) == 0
    print(io, "0")
    if finiteprecision(f)
      print(io, " + O($(string(t))^($(f.exp[end])))")
    end
  else
    for i in 1:length(f.coeff) - 1
      print(io, "$(f.coeff[i]) * $(string(t))^($(f.exp[i])) + ")
    end

    print(io, "$(f.coeff[end]) * $(string(t))^($(f.exp[length(f.coeff)]))")

    if finiteprecision(f)
      print(io, " + O($(string(t))^($(f.exp[end])))")
    end
  end
end

###############################################################################
#
#   Promotions
#
###############################################################################

Base.promote_rule{T <: Integer}(::Type{PuiseuxSeries}, ::Type{T}) = PuiseuxSeries

Base.promote_rule{T}(::Type{PuiseuxSeries{T}}, ::Type{T}) = PuiseuxSeries{T}

Base.promote_rule{T}(::Type{PuiseuxSeries{T}}, ::Type{fmpz}) = PuiseuxSeries{T}

# Newton polygon

function newton_polygon(f)
  points = Array{Tuple{fmpq, fmpq}}(0)
  for i in 0:degree(f)
     c = coeff(f, i)
     if c != 0
       push!(points, (fmpq(i), fmpq(valuation(c))))
     end
  end
  println(points)
  return lowerconvexhull(points)
end

function puiseux_expansion(f, val)
  R = parent(f)
  x = gen(R)
  t = gen(base_ring(R))
  Rc, c = PolynomialRing(R, "c")
  Rx1, x1 = PolynomialRing(Rc, "x1")
  fsubst = subst(f, t^val*(c + x1))

  fsubstzero = coeff(fsubst, 0)

  Kz, z = PolynomialRing(base_ring(base_ring(R)), "z")

  current_lowest = initial(coeff(coeff(fsubstzero, degree(fsubstzero)), 0))*z^degree(fsubstzero)
  current_v = valuation(coeff(coeff(fsubstzero, degree(fsubstzero)), 0))

  for i in degree(fsubstzero) - 1:-1:0 # c is main variable
    fsubstzerocoeff = coeff(coeff(fsubstzero, i), 0) # fubstcoeff is in (R{t})[x]

    if fsubstzerocoeff == 0
      continue
    end

    if current_v > valuation(fsubstzerocoeff)
      current_lowest = initial(fsubstzerocoeff)*z^i
      current_v = valuation(fsubstzerocoeff)
    elseif current_v == valuation(fsubstzerocoeff)
      current_lowest += initial(fsubstzerocoeff)*z^i
    end
  end
  return current_lowest
end

