###############################################################################
#
#   Error.jl : error types
#
#   Copyright (C) 2017 Johannes Schmitt
#
###############################################################################

export ErrorNotSquare, _check_is_square

type ErrorNotSquare <: Exception
   get_r::Int
   get_c::Int

   function ErrorNotSquare(gr::Int, gc::Int)
      e = new(gr, gc)
      return e
   end
   
   function ErrorNotSquare(A::MatElem)
      return ErrorSquare(rows(A), cols(A))
   end
end

function Base.showerror(io::IO, e::ErrorNotSquare)
   print(io, "Expected square matrix, but matrix is $(e.get_r) x $(e.get_c).")
end

function _check_is_square(M::MatElem)
   rows(M) != cols(M) && throw(ErrorNotSquare(rows(M), cols(M)))
   return nothing
end
