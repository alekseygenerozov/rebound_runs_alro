# import numpy as np
import subprocess

def bash_command(cmd):
    '''Run command from the bash shell'''
    process=subprocess.Popen(['/bin/bash', '-c',cmd],  stdin=subprocess.PIPE, stdout=subprocess.PIPE)
    return process.communicate()[0]


smas=[2,6]
masses=[5e-2, 5e-3]

for sma in smas:
    for mass in masses:
        s=open("config_template", "r").read()
        s=s.replace("MASS", str(mass))
        s=s.replace("SMA", str(sma))
        s=s.replace("OUT", "M{0:.1e}_a{1}".format(mass, sma))
        bash_command("mkdir -p grid/M{0:.1e}_a{1}".format(mass, sma))

        f2=open("grid/M{0:.1e}_a{1}/config".format(mass, sma), "w")
        f2.write(s)
        f2.close()

        s2=open("run_sim_aleksey.sh", "r").read()
        s2=s2.replace("OUT", "M{0:.1e}_a{1}".format(mass, sma))
        f3=open("grid/M{0:.1e}_a{1}.sh".format(mass, sma), "w")
        f3.write(s2)
        f3.close()