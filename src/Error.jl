###############################################################################
#
#   Error.jl : error types
#
#   Copyright (C) 2017 Johannes Schmitt
#
###############################################################################

export NotSquareError, _check_issquare, MatConstrError, error_dim_negative

type NotSquareError <: Exception
   get_r::Int
   get_c::Int

   function NotSquareError(gr::Int, gc::Int)
      e = new(gr, gc)
      return e
   end
   
   function NotSquareError(A::MatElem)
      return NotSquareError(rows(A), cols(A))
   end
end

function Base.showerror(io::IO, e::NotSquareError)
   print(io, "Expected square matrix, but matrix is $(e.get_r) x $(e.get_c).")
end

function _check_issquare(M::MatElem)
   rows(M) != cols(M) && throw(NotSquareError(rows(M), cols(M)))
   return nothing
end

type MatConstrError <: Exception
  expect_r::Int
  expect_c::Int
  get_r::Int
  get_c::Int
  get_l::Int

  function MatConstrError(er::Int, ec::Int, gr::Int, gc::Int)
    e = new(er, ec, gr, gc, -1)
    return e
  end

  function MatConstrError(er::Int, ec::Int, gl::Int)
    e = new(er, ec, -1, -1, gl)
    return e
  end

  function MatConstrError{T}(er::Int, ec::Int, a::Array{T, 2})
    gr, gc = size(a)
    return MatConstrError(er, ec, gr, gc)
  end

  function MatConstrError{T}(er::Int, ec::Int, a::Array{T, 1})
    gl = length(a)
    return MatConstrError(er, ec, gl)
  end
end

function Base.showerror(io::IO, e::MatConstrError)
  if e.get_l == -1
    print(io, "Expected dimension $(e.expect_r) x $(e.expect_c), ")
    print(io, "got $(e.get_r) x $(e.get_c)")
  else
    print(io, "Expected an array of length $(e.expect_r * e.expect_c), ")
    print(io, "got $(e.get_l)")
  end
end

const error_dim_negative = ErrorException("Dimensions must be non-negative")
