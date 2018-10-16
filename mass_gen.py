# from bash_command import bash_command

#bash_commmand.bash_command('ls archive*bin >names')
names=open('names', 'r')
# with open('/home/aleksey/Dropbox/projects/disk_binaries/sim_data/nstar_experiments_2/names', 'r') as names:
dat=[line.strip() for line in names]
names.close()

temp=open('/home/aleksey/Dropbox/projects/disk_binaries/sim_data/template.sh', 'r')
s0=temp.read()
temp.close()
for idx,name in enumerate(dat):
    out=open('mass_{0}.sh'.format(idx), 'w')
    s=s0+'\npython /projects/alge9397/code/python/rebound_runs/mass_script.py {0}'.format(name)
    out.write(s)
out.close()