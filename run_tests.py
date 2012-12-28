#aaa
#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2011, Marine Biological Laboratory
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.
#

from unittest import TestLoader, TextTestRunner, TestSuite
import test.illumina_files_tests as ill_f
import test.db_upload_tests as db_up
import test.utils_tests as utils

if __name__ == "__main__":

    loader = TestLoader()
    suite = TestSuite((
        loader.loadTestsFromTestCase(utils.IlluminaFilesTestCase),
        loader.loadTestsFromTestCase(ill_f.IlluminaFilesTestCase),
        loader.loadTestsFromTestCase(db_up.DbUloadTestCase),
        ))

    runner = TextTestRunner(verbosity = 2)
    runner.run(suite)
