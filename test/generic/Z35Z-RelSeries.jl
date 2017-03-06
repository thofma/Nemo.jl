# To define a powerseries ring
type Z35Z <: Nemo.Ring end
type Z35Z_elem <: Nemo.RingElem
   d::Int
   parent::Z35Z
   Z35Z_elem(a::Int) = new(a%35)
   Z35Z_elem(a::fmpz) = new(a%35)
end

Nemo.elem_type(::Z35Z) = Z35Z_elem

Nemo.parent_type(::Type{Z35Z_elem}) = Z35Z

(a::Z35Z)() = a(0)

function (a::Z35Z)(b::Int)
   c = Z35Z_elem(b)
   c.parent = a
   return c
end

Nemo.parent(a::Z35Z_elem) = a.parent

Nemo.needs_parentheses(a::Z35Z_elem) = false

Nemo.iszero(a::Z35Z_elem) = a.d%35 == 0

Nemo.isone(a::Z35Z_elem) = a.d%35 == 1 || a.d%35 == -34

Base.promote_rule(a::Type{Z35Z_elem}, b::Type{Int}) = Z35Z_elem

Base.promote_rule(a::Type{Z35Z_elem}, b::Type{Z35Z_elem}) = Z35Z_elem

function (a::Z35Z)(b::Z35Z_elem)
   c = b
   c.parent = a
   return c
end

import Base.==
function ==(a::Z35Z_elem, b::Z35Z_elem)
   a1 = (a.d+35)%35
   b1 = (b.d+35)%35
   return a1 == b1
end

# For arithmetic

Nemo.zero(a::Z35Z) = a(0)

import Base.+
function +(a::Z35Z_elem, b::Z35Z_elem)
   c = Z35Z_elem(a.d+b.d)
   c.parent = a.parent
   return c
end

import Base.-
function -(a::Z35Z_elem, b::Z35Z_elem)
   c = Z35Z_elem(a.d-b.d)
   c.parent = a.parent
   return c
end

function -(a::Z35Z_elem)
   c = Z35Z_elem(-a.d)
   c.parent = a.parent
   return c
end

Nemo.show_minus_one(::Type{Z35Z_elem}) = false

import Base.*
function *(a::Z35Z_elem, b::Z35Z_elem)
   c = Z35Z_elem(a.d*b.d)
   c.parent = a.parent
   return c
end

Nemo.is_negative(a::Z35Z_elem) = a.d < 0

import Nemo.mul!
function Nemo.mul!(c::Z35Z_elem, a::Z35Z_elem, b::Z35Z_elem)
   c.d = (a.d*b.d)%35
   c.parent = a.parent
   nothing
end

import Nemo.addeq!
function Nemo.addeq!(b::Z35Z_elem, a::Z35Z_elem)
   b.d = (b.d+a.d)%35
   nothing
end

# For test_rel_series_constructors
function (a::Z35Z)(b::fmpz)
   c = Z35Z_elem(b)
   c.parent = a
   return c
end

# For test_rel_series_manipulation

function Nemo.isunit(a::Z35Z_elem)
   if a.d % 5 == 0 || a.d % 7 == 0
      return false
   end
   return true
end

# For test_rel_series_adhoc_binary_ops

Base.promote_rule(a::Type{Z35Z_elem}, b::Type{fmpz}) = Z35Z_elem

# For test_rel_series_unsafe_ops

Nemo.one(a::Z35Z) = a(1)

# For test_rel_series_inversion

function Nemo.divexact(a::Z35Z_elem, b::Z35Z_elem)
   g, binv = gcdinv(ZZ((b.d+35)%35), ZZ(35))
   if g != 1
      error("Impossible inverse in divexact")
   end
   return parent(a)(a.d * binv)
end
