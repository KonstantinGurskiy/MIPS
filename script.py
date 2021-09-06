act = 70
while act<100:
    f=open('/run/media/softmatter/Новый том/Hoomd/in.py','r')
    lines=f.readlines()
    lines[19]='Act = '+str(act)+'\n'
    #lines[9]='hoomd.init.create_lattice(unitcell=hoomd.lattice.hex(a='+str( (2/((0.3+0.1*a)*(3**(0.5)))) ** (0.5) ) + '), n=[100,50])'+'\n'
    #  lines[24]='hoomd.md.force.active(seed=24, group=all, f_lst=activity, orientation_link = False, rotation_diff='+(str(0.1+0.1*dr))+')'+'\n'
    #  lines[30]='bd.set_gamma('+"'A'"+','+str(0.1+0.1*dt)+')'+'\n'
    lines[33]='hoomd.analyze.log(filename="'+'log-output.log'+'act'+str(act)+'", quantities=['+"'potential_energy'"+'], period=20000, overwrite=True)'+'\n'
    lines[35]='dp = hoomd.deprecated.dump.xml(group=all, filename="trajectory'+'act'+str(act)+'", period=20000)'+'\n'
    f.close()
    save_changes = open('/run/media/softmatter/Новый том/Hoomd/in.py', 'w')
    save_changes.writelines(lines)
    save_changes.close()
    act=act+5
    import os
    os.system('python3 in.py --mode=gpu')
