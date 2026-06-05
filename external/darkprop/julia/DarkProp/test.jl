include("DarkProp.jl")

rn = DarkProp.RandomNumber(-1)
@show DarkProp.uniform_phi(rn)
@show DarkProp.uniform_costheta(rn)
@show DarkProp.uniform_xi(rn)
@show DarkProp.uniform_ab(rn, -10, 10.0)
@show DarkProp.normal(rn, 0, 1)

O16 = DarkProp.Target("Oxygen", 8, 16, 15.0)
empty = DarkProp.Target()
@show DarkProp.get_Z(O16)
@show DarkProp.get_A(O16)
@show DarkProp.get_name(O16)
@show DarkProp.get_m(O16)
@show DarkProp.get_A(empty)

v1 = DarkProp.Vector3d(1.0, 2.0, 3.0)
v2 = DarkProp.Vector3d(3.0, 2.0, 0)
@show vec(v1)
@show vec(v2)
@show vec(-v2)
@show vec(v1 - v2)
@show vec(v2 - v1)
@show vec(v1 + v2)
@show vec(v2 + v1)
@show v1 == v2
@show typeof(vec(DarkProp.Vector3d(1.0, 2, 3)))

dm = DarkProp.SIDM()
@show DarkProp.in_earth(dm)
DarkProp.in_earth!(dm, true)
@show DarkProp.in_earth(dm)
DarkProp.in_earth!(dm, false)
@show DarkProp.in_earth(dm)

@show dm isa DarkProp.Particle

earth = DarkProp.HomoEarth()
@show earth

DarkProp.scatter!(dm, earth, rn)
