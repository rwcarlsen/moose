#!/usr/bin/env python
"""
python /Users/calsrw/moose/examples/ex04_bcs/dirichlet_convected.py
"""

import vtk
import chigger

def makechart(infile, outfile, varlabel, varname, lim, t=0):
    camera = vtk.vtkCamera()
    camera.SetViewUp(0.7247, -0.2244, 0.6515)
    camera.SetPosition(-5.9077, 4.2951, 7.3355)
    camera.SetFocalPoint(-0.5513, -0.2652, -0.1935)

    reader = chigger.exodus.ExodusReader(infile, timestep=None, time=t)
    reader.setOptions(block=['1'])

    result = chigger.exodus.ExodusResult(reader)
    result.setOptions(edges=True, edge_color=[0, 0, 0], variable=varname, block=['1'], min=lim[0], max=lim[1], cmap='plasma', local_range=True, camera=camera)

    cbar = chigger.exodus.ExodusColorBar(result)
    cbar.setOptions(colorbar_origin=(0.8, 0.25, 0.0), title='{} ({})'.format(varlabel, varname))
    cbar.setOptions('primary', lim=lim, font_size=20, font_color=[0, 0, 0])

    #extents = chigger.misc.VolumeAxes(result)
    #extents.setOptions('xaxis', 'yaxis', color=[0,0,0])

    window = chigger.RenderWindow(result, cbar)
    window.setOptions(background=[1, 1, 1], size=[1000, 1000], style='interactive', antialiasing=200, test=True)
    window.write(outfile)

if __name__ == '__main__':
    makechart('transient_out.e', 'transient_t0.png', 'u', 'diffused', [0, 8], t=0)
    makechart('transient_out.e', 'transient_tmid.png', 'u', 'diffused', [0, 8], t=.1)
    makechart('transient_out.e', 'transient_tend.png', 'u', 'diffused', [0, 8], t=1)
