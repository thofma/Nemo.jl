include("Z35Z-Poly.jl")

function test_gen_poly_constructors()
   print("GenPoly.constructors...")
 
   R, x = ZZ["x"]
   S, y = R["y"]

   @test typeof(R) <: Nemo.Ring
   @test typeof(S) <: GenPolyRing

   @test isa(y, PolyElem)

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test typeof(S) <: GenPolyRing

   @test isa(y, PolyElem)

   T, z = PolynomialRing(S, "z")

   @test typeof(T) <: GenPolyRing

   @test isa(z, PolyElem)

   f = x^2 + y^3 + z + 1

   @test isa(f, PolyElem)

   g = S(2)

   @test isa(g, PolyElem)

   h = S(x^2 + 2x + 1)

   @test isa(h, PolyElem)

   j = T(x + 2)

   @test isa(j, PolyElem)

   k = S([x, x + 2, x^2 + 3x + 1])

   @test isa(k, PolyElem)

   l = S(k)

   @test isa(l, PolyElem)

   U, u = PolynomialRing(Z35Z(), "u")

   @test typeof(U) <: GenPolyRing

   @test isa(u, PolyElem)

   f = u^2 + 1

   @test isa(f, PolyElem)

   g = U(2)

   @test isa(f, PolyElem)

   println("PASS")
end

function test_gen_poly_manipulation()
   print("GenPoly.manipulation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test iszero(zero(S))
   
   @test isone(one(S))

   @test isgen(gen(S))
   
   @test isunit(one(S))

   f = 2x*y + x^2 + 1

   @test lead(f) == 2x

   @test degree(f) == 1

   h = x*y^2 + (x + 1)*y + 3

   @test coeff(h, 2) == x

   @test length(h) == 3

   @test canonical_unit(-x*y + x + 1) == -1

   @test deepcopy(h) == h

   U, u = PolynomialRing(Z35Z(), "u")

   @test iszero(zero(U))
   
   @test isone(one(U))

   @test isgen(gen(U))
   
   @test isunit(one(U))

   f = 4*u^2 + 1

   @test lead(f) == U(4)

   @test degree(f) == 2

   h = 2*u^2 + 4*u + 3

   @test coeff(h, 2) == 2

   @test length(h) == 3

   @test canonical_unit(-u^2+1) == -1

   @test deepcopy(h) == h

   println("PASS")
end

function test_gen_poly_binary_ops()
   print("GenPoly.binary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)
   
   @test f - g == x*y^2+(-x^3-2*x+1)

   @test f + g == x*y^2+(2*x+2)*y+(x^3+2*x+5)

   @test f*g == (x^2+x)*y^3+(x^4+3*x^2+4*x+1)*y^2+(x^4+x^3+2*x^2+7*x+5)*y+(3*x^3+6*x+6)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2
   
   @test f - g == -2u+1

   @test f + g == 2u^3+2u^2+4u+5

   @test f*g == u^6+2u^5+5u^4+9u^3+8u^2+11u+6 

   println("PASS")
end

function test_gen_poly_adhoc_binary()
   print("GenPoly.adhoc_binary...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)

   @test f*4 == (4*x)*y^2+(4*x+4)*y+12

   @test 7*f == (7*x)*y^2+(7*x+7)*y+21
   
   @test fmpz(5)*g == (5*x+5)*y+(5*x^3+10*x+10)

   @test g*fmpz(3) == (3*x+3)*y+(3*x^3+6*x+6)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test f*4 == 4u^3+4u^2+4u+12

   @test 7*f == 7u^3+7u^2+7u+21
   
   @test fmpz(5)*g == 5u^3+5u^2+15u+10

   @test g*fmpz(3) == 3u^3+3u^2+9u+6

   println("PASS")
end

function test_gen_poly_unsafe_ops()
   print("GenPoly.unsafe_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   zero!(f)
   @test iszero(f)

   one!(f)
   @test isone(f)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3+4

   zero!(f)
   @test iszero(f)

   one!(f)
   @test isone(f)

   println("PASS")
end

function test_gen_poly_comparison()
   print("GenPoly.comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = x*y^2 + (x + 1)*y + 3

   @test f == g

   @test isequal(f, g)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3+u^2+3
   g = u^3+u^2+3

   @test f == g

   @test isequal(f, g)

   println("PASS")
end

function test_gen_poly_adhoc_comparison()
   print("GenPoly.adhoc_comparison...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test S(1) == 1 

   @test 1 != x + y

   U, u = PolynomialRing(Z35Z(), "u")

   @test U(1) == 1 

   @test 1 != 2u

   println("PASS")
end

function test_gen_poly_unary_ops()
   print("GenPoly.unary_ops...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test -f == -x*y^2 - (x + 1)*y - 3

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + 3u - 1

   @test -f == -u^3-3u+1

   println("PASS")
end

function test_gen_poly_truncation()
   print("GenPoly.truncation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)

   @test truncate(f, 1) == 3

   @test mullow(f, g, 4) == (x^2+x)*y^3+(x^4+3*x^2+4*x+1)*y^2+(x^4+x^3+2*x^2+7*x+5)*y+(3*x^3+6*x+6)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test truncate(f, 1) == 3

   @test mullow(f, g, 4) == 9u^3+8u^2+11u+6

   println("PASS")
end

function test_gen_poly_reverse()
   print("GenPoly.reverse...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test reverse(f, 7) == 3y^6 + (x + 1)*y^5 + x*y^4

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + 2u^2 + u + 3

   @test reverse(f, 7) == 3u^6 + u^5 + 2u^4 + u^3

   println("PASS")
end

function test_gen_poly_shift()
   print("GenPoly.shift...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test shift_left(f, 7) == x*y^9 + (x + 1)*y^8 + 3y^7

   @test shift_right(f, 3) == 0

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + 2u^2 + u + 3

   @test shift_left(f, 7) == u^10 + 2u^9 + u^8 + 3u^7

   @test shift_right(f, 4) == 0

   println("PASS")
end

function test_gen_poly_powering()
   print("GenPoly.powering...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test f^5 == (x^5)*y^10+(5*x^5+5*x^4)*y^9+(10*x^5+35*x^4+10*x^3)*y^8+(10*x^5+90*x^4+90*x^3+10*x^2)*y^7+(5*x^5+110*x^4+300*x^3+110*x^2+5*x)*y^6+(x^5+65*x^4+460*x^3+460*x^2+65*x+1)*y^5+(15*x^4+330*x^3+900*x^2+330*x+15)*y^4+(90*x^3+810*x^2+810*x+90)*y^3+(270*x^2+945*x+270)*y^2+(405*x+405)*y+243

   g = shift_left(f, 3)

   @test g^5 == y^15*f^5

   U, u = PolynomialRing(Z35Z(), "u")

   f = 2u^2 + u + 3

   @test f^5 == 32*u^10+10*u^9+5*u^8+30*u^7+5*u^6+11*u^5+25*u^4+15*u^3+30*u^2+20*u+33 

   g = shift_left(f, 3)

   @test g^5 == u^15*f^5

   println("PASS")
end

function test_gen_poly_modular_arithmetic()
   print("GenPoly.modular_arithmetic...")

   R, x = PolynomialRing(QQ, "x")
   S = ResidueRing(R, x^3 + 3x + 1)
   T, y = PolynomialRing(S, "y")

   f = (3*x^2 + x + 2)*y + x^2 + 1
   g = (5*x^2 + 2*x + 1)*y^2 + 2x*y + x + 1
   h = (3*x^3 + 2*x^2 + x + 7)*y^5 + 2x*y + 1

   @test invmod(f, g) == (fmpz(707)//3530*x^2 + fmpz(2151)//1765*x + fmpz(123)//3530)*y+(fmpz(-178)//1765*x^2 - fmpz(551)//3530*x + fmpz(698)//1765)
   
   @test mulmod(f, g, h) == (-30*x^2 - 43*x - 9)*y^3+(-7*x^2 - 23*x - 7)*y^2+(4*x^2 - 10*x - 3)*y+(x^2 - 2*x)
   
   @test powmod(f, 3, h) == (69*x^2 + 243*x + 79)*y^3+(78*x^2 + 180*x + 63)*y^2+(27*x^2 + 42*x + 18)*y+(3*x^2 + 3*x + 2)

   println("PASS")
end

function test_gen_poly_exact_division()
   print("GenPoly.exact_division...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)

   @test divexact(f*g, f) == g

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test divexact(f*g, f) == g

   println("PASS")
end

function test_gen_poly_adhoc_exact_division()
   print("GenPoly.adhoc_exact_division...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test divexact(3*f, 3) == f

   @test divexact(x*f, x) == f

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3

   @test divexact(3*f, 3) == f

   @test divexact(u*f, u) == f

   println("PASS")
end

function test_gen_poly_euclidean_division()
   print("GenPoly.euclidean_division...")

   R = ResidueRing(ZZ, 7)
   S, x = PolynomialRing(R, "x")
   T = ResidueRing(S, x^3 + 3x + 1)
   U, y = PolynomialRing(T, "y")

   k = y^3 + x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y^2 + (x^3 + 2x + 2)

   @test mod(k, l) == (4*x^2+4*x+4)*y+(3*x^2+5*x+6)
   
   @test divrem(k, l) == ((5*x^2+2*x+6)*y+(2*x^2+5*x+2), (4*x^2+4*x+4)*y+(3*x^2+5*x+6))
 
   println("PASS")
end

function test_gen_poly_pseudodivision()
   print("GenPoly.pseudodivision...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)

   @test pseudorem(k, l) == (x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1)

   @test pseudodivrem(k, l) == ((x^2+x)*y+(-x^4-x^2+1), (x^7+3*x^5+2*x^4+x^3+5*x^2+4*x+1))

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test pseudorem(f, g) == 33u+1

   @test pseudodivrem(f, g) == (1, 33u+1)

   println("PASS")
end

function test_gen_poly_content_primpart_gcd()
   print("GenPoly.content_primpart_gcd...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   k = x*y^2 + (x + 1)*y + 3
   l = (x + 1)*y + (x^3 + 2x + 2)
   m = y^2 + x + 1

   @test content(k) == 1

   @test primpart(k*(x^2 + 1)) == k

   @test gcd(k*m, l*m) == m

   @test lcm(k*m, l*m) == k*l*m

   R, x = PolynomialRing(QQ, "x")
   T = ResidueRing(R, x^3 + 3x + 1)
   U, z = PolynomialRing(T, "z")

   r = z^3 + 2z + 1
   s = z^5 + 1

   @test gcdinv(r, s) == (1,fmpz(-21)//62*z^4+fmpz(13)//62*z^3-fmpz(11)//62*z^2-fmpz(5)//62*z+fmpz(9)//62)

   U, u = PolynomialRing(Z35Z(), "u")

   k = u^3 + u^2 + u + 3

   @test content(k) == 1

   @test primpart(k*(u^2 + 1)) == u^5+u^4+2*u^3+4*u^2+u+3

   println("PASS")
end

function test_gen_poly_evaluation()
   print("GenPoly.evaluation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x^2 + 2x + 1
   g = x*y^2 + (x + 1)*y + 3

   @test evaluate(g, 3) == 12x + 6

   @test evaluate(g, fmpz(3)) == 12x + 6

   @test evaluate(g, f) == x^5+4*x^4+7*x^3+7*x^2+4*x+4

   U, u = PolynomialRing(Z35Z(), "u")

   g = u^3 + u^2 + 3u + 2

   @test evaluate(g, 3) == 12

   @test evaluate(g, fmpz(3)) == 12

   @test evaluate(g, Z35Z()(3)) == 12 

   println("PASS")
end

function test_gen_poly_composition()
   print("GenPoly.composition...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3
   g = (x + 1)*y + (x^3 + 2x + 2)

   @test compose(f, g) == (x^3+2*x^2+x)*y^2+(2*x^5+2*x^4+4*x^3+9*x^2+6*x+1)*y+(x^7+4*x^5+5*x^4+5*x^3+10*x^2+8*x+5)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test compose(f, g) == u^9+3u^8+12u^7+26u^6+15u^5+6u^4+16u^3+10u^2+16u+17

   println("PASS")
end

function test_gen_poly_derivative()
   print("GenPoly.derivative...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   h = x*y^2 + (x + 1)*y + 3

   @test derivative(h) == 2x*y + x + 1

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3

   @test derivative(f) == 3u^2+2u+1

   println("PASS")
end

function test_gen_poly_integral()
   print("GenPoly.integral...")

   R, x = PolynomialRing(QQ, "x")
   S = ResidueRing(R, x^3 + 3x + 1)
   T, y = PolynomialRing(S, "y")

   f = (x^2 + 2x + 1)*y^2 + (x + 1)*y - 2x + 4
   
   @test integral(f) == (fmpz(1)//3*x^2 + fmpz(2)//3*x + fmpz(1)//3)*y^3+(fmpz(1)//2*x+fmpz(1)//2)*y^2+(-2*x+4)*y

   println("PASS")
end

function test_gen_poly_resultant()
   print("GenPoly.resultant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test resultant(f, g) == 3*x^7+6*x^5-6*x^3+96*x^2+192*x+96

   println("PASS")
end

function test_gen_poly_discriminant()
   print("GenPoly.discriminant...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = x*y^2 + (x + 1)*y + 3

   @test discriminant(f) == x^2-10*x+1

   println("PASS")
end

function test_gen_poly_gcdx()
   print("GenPoly.gcdx...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test gcdx(f, g) == (3*x^7+6*x^5-6*x^3+96*x^2+192*x+96, (36*x^2+72*x+36), (-18*x^2-18*x)*y+(3*x^4-6*x-6))

   println("PASS")
end

function test_gen_poly_newton_representation()
   print("GenPoly.newton_representation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3

   g = deepcopy(f)
   roots = [R(1), R(2), R(3)]
   monomial_to_newton!(g.coeffs, roots)
   newton_to_monomial!(g.coeffs, roots)

   @test f == g

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3

   g = deepcopy(f)
   roots = [Z35Z()(1), Z35Z()(2), Z35Z()(3)]
   monomial_to_newton!(g.coeffs, roots)
   newton_to_monomial!(g.coeffs, roots)

   @test f == g

   println("PASS")
end

function test_gen_poly_interpolation()
   print("GenPoly.interpolation...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   xs = [R(1), R(2), R(3), R(4)]
   ys = [R(1), R(4), R(9), R(16)]

   f = interpolate(S, xs, ys)

   @test f == y^2

   U, u = PolynomialRing(Z35Z(), "u")

   xs = [Z35Z()(1), Z35Z()(2), Z35Z()(3), Z35Z()(4)]
   ys = [Z35Z()(1), Z35Z()(4), Z35Z()(9), Z35Z()(16)]

   f = interpolate(U, xs, ys)

   @test f == u^2

   println("PASS")
end

function test_gen_poly_special()
   print("GenPoly.special...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   @test chebyshev_t(20, y) == 524288*y^20-2621440*y^18+5570560*y^16-6553600*y^14+4659200*y^12-2050048*y^10+549120*y^8-84480*y^6+6600*y^4-200*y^2+1

   @test chebyshev_u(15, y) == 32768*y^15-114688*y^13+159744*y^11-112640*y^9+42240*y^7-8064*y^5+672*y^3-16*y

   U, u = PolynomialRing(Z35Z(), "u")

   @test chebyshev_t(20, u) == 23u^20+25u^18+30u^16+10u^14+7u^10+5u^8+10u^6+20u^4+10u^2+1

   @test chebyshev_u(15, u) == 8u^15+7u^13+4u^11+25u^9+30u^7+21u^5+7u^3+19u

   println("PASS")
end

function test_gen_poly_mul_karatsuba()
   print("GenPoly.mul_karatsuba...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   
   f = x + y + 2z^2 + 1
   
   @test mul_karatsuba(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_karatsuba(f^10, f^30) == mul_classical(f^10, f^30)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3

   @test mul_karatsuba(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_karatsuba(f^10, f^30) == mul_classical(f^10, f^30)

   println("PASS")
end

function test_gen_poly_mul_ks()
   print("GenPoly.mul_ks...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")
   T, z = PolynomialRing(S, "z")
   
   f = x + y + 2z^2 + 1
   
   @test mul_ks(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_ks(f^10, f^30) == mul_classical(f^10, f^30)

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3

   @test mul_karatsuba(f^10, f^10) == mul_classical(f^10, f^10)
   @test mul_karatsuba(f^10, f^30) == mul_classical(f^10, f^30)

   println("PASS")
end

function test_gen_poly_generic_eval()
   print("GenPoly.generic_eval...")

   R, x = PolynomialRing(ZZ, "x")
   S, y = PolynomialRing(R, "y")

   f = 3x*y^2 + (x + 1)*y + 3
   g = 6(x + 1)*y + (x^3 + 2x + 2)

   @test f(g) == (108*x^3+216*x^2+108*x)*y^2+(36*x^5+36*x^4+72*x^3+150*x^2+84*x+6)*y+(3*x^7+12*x^5+13*x^4+13*x^3+26*x^2+16*x+5)

   @test f(x + 1) == 3*x^3+7*x^2+5*x+4

   @test f(123) == 45510*x+126

   @test f(fmpz(123)) == 45510*x + 126

   R, x = PolynomialRing(ZZ, "x")
   T, y = FiniteField(103, 1, "y")

   f = x^5 + 3x^3 + 2x^2 + x + 1

   @test f(T(13)) == 20

   U, u = PolynomialRing(Z35Z(), "u")

   f = u^3 + u^2 + u + 3
   g = u^3 + u^2 + 3u + 2

   @test f(g) == u^9+3u^8+12u^7+26u^6+15u^5+6u^4+16u^3+10u^2+16u+17

   @test f(u + 1) == u^3+4u^2+6u+6 

   @test f(123) == 17

   @test f(fmpz(123)) == 17

   println("PASS")
end

function test_gen_poly()
   test_gen_poly_constructors()
   test_gen_poly_manipulation()
   test_gen_poly_binary_ops()
   test_gen_poly_adhoc_binary()
   test_gen_poly_unsafe_ops()
   test_gen_poly_comparison()
   test_gen_poly_adhoc_comparison()
   test_gen_poly_unary_ops()
   test_gen_poly_truncation()
   test_gen_poly_reverse()
   test_gen_poly_shift()
   test_gen_poly_powering()
   test_gen_poly_modular_arithmetic()
   test_gen_poly_exact_division()
   test_gen_poly_adhoc_exact_division()
   test_gen_poly_euclidean_division()
   test_gen_poly_pseudodivision()
   test_gen_poly_content_primpart_gcd()
   test_gen_poly_evaluation()
   test_gen_poly_composition()
   test_gen_poly_derivative()
   test_gen_poly_integral()
   test_gen_poly_resultant()
   test_gen_poly_discriminant()
   test_gen_poly_gcdx()
   test_gen_poly_newton_representation()
   test_gen_poly_interpolation()
   test_gen_poly_special()
   test_gen_poly_mul_karatsuba()
   test_gen_poly_mul_ks()
if VERSION >= v"0.5.0-dev+3171"
   test_gen_poly_generic_eval()
end

   println("")
end
