#!/usr/bin/env python3

"""writeXDMF.py:
Script to parse an HDF5 file written by postprocessing Sullivan's ABL data and write
corresponding XDMF file(s), allowing the data to be postprocessed.
"""



import h5py
import numpy as np
import math
import re
import argparse
import os

# Parse and write mesh/fields
def writeFields(fo):
    
    # Print header
    fo.write('<Xdmf>\n')
    fo.write('  <Domain>\n')
    fo.write('    <Grid Name="FieldData" GridType="Collection" '
             'CollectionType="Temporal">\n\n')

    # Loop over all times
    for index in range(np.size(timeDirs)):

        timeName = timeDirs[index]


        fo.write('      <Grid Name="mesh1" GridType="Uniform">\n')
        fo.write('        <Time Type="Single" Value="{}" />\n'.format(timeName-38219.4))
        fo.write('        <Topology TopologyType="3DCoRectMesh" NumberOfElements="50 768 768"/>\n')
        fo.write('        <Geometry GeometryType="ORIGIN_DXDYDZ">\n')
        fo.write('        <DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n')
        fo.write('        0 0 0\n')
        fo.write('        </DataItem>\n')
        fo.write('        <DataItem Dimensions="3" NumberType="Float" Precision="4" Format="XML">\n')
        fo.write('        8.0 6.6667 6.6667\n')
        fo.write('        </DataItem>\n')
        fo.write('        </Geometry>\n')
        fo.write('        <Attribute Name="U" AttributeType="Vector" Center="Cell">\n')
        fo.write('        <DataItem ItemType="Function"  Function="JOIN($0,$1,$2)" Dimensions="50 768 768 3">\n')
        fo.write('        <DataItem Dimensions="50 768 768" NumberType="Float" Precision="4" Format="HDF">\n')
        fo.write('        velData_{}.h5:/Ux\n'.format(timeName))
        fo.write('        </DataItem>\n')
        fo.write('        <DataItem Dimensions="50 768 768" NumberType="Float" Precision="4" Format="HDF">\n')
        fo.write('        velData_{}.h5:/Uy\n'.format(timeName))
        fo.write('        </DataItem>\n')
        fo.write('        <DataItem Dimensions="50 768 768" NumberType="Float" Precision="4" Format="HDF">\n')
        fo.write('        velData_{}.h5:/Uz\n'.format(timeName))
        fo.write('        </DataItem>\n')
        fo.write('        </DataItem>\n')
        fo.write('        </Attribute>\n')
        fo.write('        </Grid>     \n')

    # Print footer
    fo.write('    </Grid>\n')
    fo.write('  </Domain>\n')
    fo.write('</Xdmf>\n');
    
    # Return control
    return



""" Main part of program/script start below here
"""
# Set up and read command line arguments
parser = argparse.ArgumentParser(description='Script to parse an HDF5 file '
                                 'written by OpenFOAM  and write corresponding '
                                 'XDMF files')
parser.add_argument('-l','--noLagrangian',
                    action='store_true',
                    help='Skip Lagrangian clouds',
                    required=False)
parser.add_argument('-f','--noFields',
                    action='store_true',
                    help='Skip fields (mesh) data',
                    required=False)
parser.add_argument('-z','--noZero',
                    action='store_true',
                    help='Do not process time zero (t=0)',
                    required=False)
parser.add_argument('-t','--latestTime',
                    action='store_true',
                    help='Only process latest timestep present',
                        required=False)
# parser.add_argument('--time',
#                     help='Process times from t1:t2',
#                         required=True)

args = parser.parse_args()

print "Reading time list from 'timeFiles'. HDF5 file names are assumed to be velData_timeName.h5"
timeDirs = np.loadtxt('timeFiles')

with open('ablData.xdmf','w') as fo:
    writeFields(fo)

