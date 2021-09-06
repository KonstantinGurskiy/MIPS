import hoomd
import hoomd.md
import hoomd.deprecated
import math
import hoomd.dem
import numpy

# initialize
hoomd.context.initialize()
hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a=1.2), n=[100,50])

# specify potential
nl = hoomd.md.nlist.cell()
lj = hoomd.md.pair.lj(r_cut=2.5, nlist=nl)
lj.pair_coeff.set('A', 'A', epsilon=0.6, sigma=1.0)

#specify force
all = hoomd.group.all();
N = len(all);
Act = 90
theta = [2.0*math.pi*numpy.random.random_sample() for i in range(N)]
activity = [(Act*math.cos(1), Act*math.sin(1),0) for i in range(N)]
for i in range(N):
    activity[i] = (Act*math.cos(theta[i]), Act*math.sin(theta[i]),0)
hoomd.md.force.active(seed=24, group=all, f_lst=activity, orientation_link = False, rotation_diff=0.05)

# define integrator
all = hoomd.group.all();
hoomd.md.integrate.mode_standard(dt=0.001)
bd = hoomd.md.integrate.langevin(group=all, kT=0.1, seed=42)
bd.set_gamma('A',2)

# write output
hoomd.analyze.log(filename="log-output.logact95", quantities=['potential_energy'], period=20000, overwrite=True)
#hoomd.dump.gsd("trajectory.gsd", period=2e3, group=all, overwrite=True)
dp = hoomd.deprecated.dump.xml(group=all, filename="trajectoryact95", period=20000)
dp.set_params(velocity=True)

# run simulation
hoomd.run(1000000)
