# -*- python -*-
# $RCSfile: cleanup.py,v $
# $Revision: 1.2 $
# $Author: langer $
# $Date: 2014/09/27 21:42:13 $

# This software was produced by NIST, an agency of the U.S. government,
# and by statute is not subject to copyright in the United States.
# Recipients of this software assume all responsibilities associated
# with its operation, modification and maintenance. However, to
# facilitate maintenance we ask that before distributing modified
# versions of this software, you first contact the authors at
# oof_manager@nist.gov. 

# This file is run with execfile from within guitests.py.  "dir" is
# set to the directory containing the test files.

import filecmp
assert filecmp.cmp("quit.log", os.path.join(dir, "quit.log"))
removefile("quit.log")
