#!/usr/bin/env python
"""
python /Users/calsrw/moose/examples/ex04_bcs/dirichlet_convected.py
"""

import vtk
import chigger

def makechart(infile, outfile, varlabel, varname, lim):
    camera = vtk.vtkCamera()
    camera.SetViewUp(0.0000, 1.0000, 0.0000)
    camera.SetPosition(0.5978, 0.4574, 2.1436)
    camera.SetFocalPoint(0.5978, 0.4574, -0.0000)

    reader = chigger.exodus.ExodusReader(infile)
    reader.setOptions(block=['1'])

    result = chigger.exodus.ExodusResult(reader)
    result.setOptions(edges=True, edge_color=[0, 0, 0], variable=varname, block=['1'], min=lim[0], max=lim[1], cmap='plasma', local_range=True, camera=camera)

    cbar = chigger.exodus.ExodusColorBar(result)
    cbar.setOptions(colorbar_origin=(0.8, 0.25, 0.0), title='{} ({})'.format(varlabel, varname))
    cbar.setOptions('primary', lim=lim, font_size=20, font_color=[0, 0, 0])

    extents = chigger.misc.VolumeAxes(result)
    extents.setOptions('xaxis', 'yaxis', color=[0,0,0])

    window = chigger.RenderWindow(result, cbar, extents)
    window.setOptions(background=[1, 1, 1], size=[1250, 1000], style='interactive', antialiasing=100, test=True)
    window.write(outfile)


import matplotlib.pyplot as plt

def makelineplot(title, infile, outfile, varlabels, varnames, ylims):
    reader = chigger.exodus.ExodusReader(infile)

    plt.clf()
    plt.cla()
    plt.title(title)
    plt.ylabel('Solution')
    plt.xlabel('Position (X)')

    lines = []
    for label, var in zip(varlabels, varnames):
        mug = chigger.exodus.ExodusResult(reader, variable=var)
        mug.update()

        p0 = (0, 0.5, 0)
        p1 = (1, 0.5, 0)
        sample = chigger.exodus.ExodusResultLineSampler(mug, point1=p0, point2=p1, resolution=200)
        sample.update()
        x = sample[0].getDistance()
        y = sample[0].getSample(var)

        plt.plot(x, y, '-', label='{} ({})'.format(var, label))
    plt.ylim(*ylims)
    plt.legend()
    plt.savefig(outfile)

if __name__ == '__main__':
    makechart('dirichlet_bc_out.e', 'dirichlet_convected.png', 'u', 'convected', [0, 2])
    makechart('dirichlet_bc_out.e', 'dirichlet_diffused.png', 'v', 'diffused', [0, 2])
    makechart('neumann_bc_out.e', 'neumann_convected.png', 'u', 'convected', [0, 1.35])
    makechart('neumann_bc_out.e', 'neumann_diffused.png', 'v', 'diffused', [0, 1.35])

    makelineplot('Dirichlet Cross-Section', 'dirichlet_bc_out.e', 'dirichlet_xsec.svg', ['u', 'v'], ['convected', 'diffused'], [0, 2])
    makelineplot('Neumann Cross-Section', 'neumann_bc_out.e', 'neumann_xsec.svg', ['u', 'v'], ['convected', 'diffused'], [0, 1.35])
