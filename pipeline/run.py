# -*- coding: utf-8 -*-
#
# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#
import sys,os
#sys.path.append('/bioware/linux/seqinfo/bin/python_pipeline/py_mbl_sequencing_pipeline')
sys.path.append(os.path.join(os.path.dirname(__file__),'..'))
from pipeline.runconfig import RunConfig

class Run(RunConfig):
    def __init__(self, config_file_path, basepythondir):
        RunConfig.__init__(self, config_file_path, basepythondir)
