nx=272#change
ny=272
nz=272
npoints=nx*ny*nz
directory = "./data/data_256/"
output = "./vtk/"
nstep=1800 
m = 0
for istep in range(1700,nstep,100):
    # read origin dat file
    # eat1_.dat
    m = istep/100
    if (m < 10):
    	fname = directory + "f_" + "7_000000" + ".dat"
    else:
    	fname = directory + "f_" + "7_000000" + ".dat"
    fp = open(fname, "r")
    data1 = fp.readlines()
    fp.close()

    # write to vtk file
    fname = output + "time_" + str(istep) + ".vtk"
    fp = open(fname, "w")

    # header of vtk file
    fp.write('# vtk DataFile Version 2.0\n')
    fp.write('time_{0:d}.vtk\n'.format(istep))
    fp.write('ASCII\n')
    fp.write('DATASET STRUCTURED_GRID\n')

    # coords of grid points
    fp.write('DIMENSIONS {0:>5d}  {1:>5d}  {2:>5d}\n'.format(nx,ny,nz))#change
    fp.write('POINTS {0:>7d} int\n'.format(npoints))#change
    #for i in range(0,nx):
    #    for j in range(0,ny):
    #        for k in range(0,nz):
    #            fp.write('{0:>14.6e}   {1:>14.6e}   {2:>14.6e}\n'.format(1.0*i,1.0*j,1.0*k))
    for i in range(0,nx):
        for j in range(0,ny):
            for k in range(0,nz):
                fp.write('{0:>d}   {1:>d}   {2:>d}\n'.format(i,j,k))

    # write grid point values:
    fp.write('POINT_DATA {0:>5d}\n'.format(npoints))#change

    fp.write('SCALARS eta12  int  1\n')
    fp.write('LOOKUP_TABLE default\n')
    for s in data1:
        v = int(s.strip())
        fp.writelines("{0:>d}\n".format(v))

    fp.close()





