# -*- coding: utf-8 -*-
# Copyright (C) 2017 Garth N. Wells
#
# This file is part of DOLFIN.
#
# DOLFIN is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# DOLFIN is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with DOLFIN. If not, see <http://www.gnu.org/licenses/>.

import sys
import os
import shutil
import zipfile

# extensions of files which must be copied with the demo when building docs
file_extensions = [".png",".gif", ".geo", ".xml", ".msh", ".pdf", ".mfront", ".sh"]
zipfile_extensions = [".geo", ".msh", ".xml", ".py", ".ipynb", ".mfront", ".sh"]
# files which must be added to zip folder
source_files_to_zip = ["shell.xdmf", "shell.h5"]

def process():
    """Copy demo rst files (C++ and Python) from the DOLFIN source tree
    into the demo source tree, and process file with pylit.
    """

    # Check that we can find pylint.py for converting foo.py.rst to
    # foo.py
    pylit_parser = "../utils/pylit/pylit.py"
    if not os.path.isfile(pylit_parser):
        raise RuntimeError("Cannot find pylit.py")

    # Directories to scan
    subdirs = ["../examples"]
    retval = os.getcwd()
    # Iterate over subdirectories containing demos
    for subdir in subdirs:

        # Get list of demos (demo name , subdirectory)
        demos = [(dI, os.path.join(subdir, dI)) for dI in os.listdir(subdir) if os.path.isdir(os.path.join(subdir, dI))]

        # Iterate over demos
        for demo, path in demos:

            # Process C++ and Python versions
            for version in ("."):

                # Get path to demo source directory (cpp/python) and
                # check that it exists
                version_path = os.path.join(path, version)
                if not os.path.isdir(version_path):
                    continue

                # Create directory in documentation tree for demo
                demo_dir = os.path.join('./demo/', demo)
                if not os.path.exists(demo_dir):
                    os.makedirs(demo_dir)

                ipynb_files =  [f for f in os.listdir(version_path) if f.endswith("ipynb") ]
                for f in ipynb_files:
                    ff = os.path.join(path, f)
                    shutil.copy(ff, demo_dir)
                    ret = os.system("jupyter-nbconvert --to script %s --output %s" % (ff,  os.path.join("../../doc",demo_dir, os.path.splitext(f)[0])))
                    if not ret == 0:
                        raise RuntimeError("Unable to convert ipynb file to a .py ({})".format(f))

                # Build list of rst files in demo source directory
                rst_files = [f for f in os.listdir(version_path) if f.endswith("rst") ]

                # Copy rst files into documentation demo directory and process with Pylit
                for f in rst_files:
                    shutil.copy(os.path.join(version_path, f), demo_dir)

                    # Run pylit on cpp.rst and py.rst files (file with 'double extensions')
                    if "py" in f.split("."):
                        rst_file = os.path.join(demo_dir, f)
                        command = pylit_parser + " " + rst_file
                        ret = os.system("python "+command)
                        if not ret == 0:
                            raise RuntimeError("Unable to convert rst file to a .py ({})".format(f))



                tocopy_files = [f for f in os.listdir(version_path) if (os.path.splitext(f)[1] in file_extensions) or f in source_files_to_zip]
                tocopy_files += [f for f in os.listdir(version_path)]
                for f in tocopy_files:
                    source = os.path.join(version_path, f)
                    print("Copying {} to {}".format(source, demo_dir))
                    shutil.copy(source, demo_dir)

                for f in ipynb_files + rst_files:
                    cwd = os.getcwd()
                    os.chdir(demo_dir)
                    files = os.listdir(".")
                    zip_files = []
                    for ff in files:
                        if any([ff.endswith(ext) for ext in zipfile_extensions]) or ff in source_files_to_zip:
                            zip_files.append(ff)
                    with zipfile.ZipFile(f.split(".")[0]+".zip", 'w') as myzip:
                        for a in zip_files:
                            myzip.write(a, compress_type=zipfile.ZIP_DEFLATED)
                    os.chdir(cwd)
