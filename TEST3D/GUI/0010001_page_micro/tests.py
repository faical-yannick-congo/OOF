# -*- python -*-
# $RCSfile: tests.py,v $
# $Revision: 1.1.2.1 $
# $Author: fyc $
# $Date: 2013/07/23 14:34:54 $

# This software was produced by NIST, an agency of the U.S. government,
# and by statute is not subject to copyright in the United States.
# Recipients of this software assume all responsibilities associated
# with its operation, modification and maintenance. However, to
# facilitate maintenance we ask that before distributing modified
# versions of this software, you first contact the authors at
# oof_manager@nist.gov. 

from generics import *

def sensitization1():
    return sensitizationCheck(
        {'Microstructure' : 1,
         'New' : 1,
         'NewFromFile' : 1,
         'Rename' : 1,
         'Copy' : 1,
         'Delete' : 1,
         'Save' : 1,
         'Pane:VoxelGroups:New' : 1,
         'Pane:VoxelGroups:Rename' : 0,
         'Pane:VoxelGroups:Copy' : 0,
         'Pane:VoxelGroups:Delete' : 0,
         'Pane:VoxelGroups:Meshable' : 0,
         'Pane:VoxelGroups:Add' : 0,
         'Pane:VoxelGroups:Remove' : 0,
         'Pane:VoxelGroups:Clear' : 0,
         'Pane:VoxelGroups:Info' : 0
         },
        base='OOF3D:Microstructure Page')

def meshableButtonState():
    return gtklogger.findWidget(
        'OOF3D:Microstructure Page:Pane:VoxelGroups:Meshable').get_active()

def sensitization2():
    return sensitizationCheck(
        {'Microstructure' : 1,
         'New' : 1,
         'NewFromFile' : 1,
         'Rename' : 1,
         'Copy' : 1,
         'Delete' : 1,
         'Save' : 1,
         'Pane:VoxelGroups:New' : 1,
         'Pane:VoxelGroups:Rename' : 1,
         'Pane:VoxelGroups:Copy' : 1,
         'Pane:VoxelGroups:Delete' : 1,
         'Pane:VoxelGroups:Meshable' : 1,
         'Pane:VoxelGroups:Add' : 1,
         'Pane:VoxelGroups:Remove' : 1,
         'Pane:VoxelGroups:Clear' : 0,
         'Pane:VoxelGroups:Info' : 1
         },
        base='OOF3D:Microstructure Page')

def sensitization3():
    return sensitizationCheck(
        {'Microstructure' : 1,
         'New' : 1,
         'NewFromFile' : 1,
         'Rename' : 1,
         'Copy' : 1,
         'Delete' : 1,
         'Save' : 1,
         'Pane:VoxelGroups:New' : 1,
         'Pane:VoxelGroups:Rename' : 1,
         'Pane:VoxelGroups:Copy' : 1,
         'Pane:VoxelGroups:Delete' : 1,
         'Pane:VoxelGroups:Meshable' : 1,
         'Pane:VoxelGroups:Add' : 1,
         'Pane:VoxelGroups:Remove' : 1,
         'Pane:VoxelGroups:Clear' : 1,
         'Pane:VoxelGroups:Info' : 1
         },
        base='OOF3D:Microstructure Page')
