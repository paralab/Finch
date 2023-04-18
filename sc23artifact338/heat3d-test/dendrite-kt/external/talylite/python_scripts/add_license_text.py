#!/usr/bin/python
"""
This script adds the license text header to the source files in the library.
A source file is defined as a file ending with an extension in the list
named source_extensions.

The text to add is given below. After the text, a separator line is added to
mark the end of the license text. This allows the script to identify existing
license text and update it with the license text below if there is a change
needed or new files are added.
Files that already have the text below will not be changed.

To use this script, simply run in the library's base directory.
It is ok to run this on files that have already had a license added by this
script in the past.
"""

import os

# this is the text that will be inserted
license_text = '''/*
  Copyright 2014-2016 Baskar Ganapathysubramanian

  This file is part of TALYFem.

  TALYFem is free software: you can redistribute it and/or modify
  it under the terms of the GNU Lesser General Public License as
  published by the Free Software Foundation, either version 2.1 of the
  License, or (at your option) any later version.

  TALYFem is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public
  License along with TALYFem.  If not, see <http://www.gnu.org/licenses/>.
*/'''
license_lines = [line + "\n" for line in license_text.split("\n")]

# this line is placed after the license text. It is used to allow the script
# to determine which files have a license added already so that that script
# will not re-add the license to files that already have it and so that it can
# update the license in existing files if needed.
sep = "// --- end license text --- //\n"

# list of directories whose files will not have the text added
exclude_dir = ["ExternalLibs", ]

# extensions of files to add the text to
source_extensions = [".cpp", ".h"]

# first grab the list of files
file_list = []

# walk through directory tree, starting in the current directory
for root, dirs, files in os.walk("."):
  # skip excluded directories
  dirs[:] = [d for d in dirs if d not in exclude_dir]
  # loop over files in directory
  for filename in files:
    add_file = False
    # check whether this file has a relevant source extension
    for extension in source_extensions:
      if filename.endswith(extension):
        add_file = True
    # add file to the list of files to modify
    if add_file:
      file_list.append(os.path.join(root, filename))

# add text to files we found
for filename in file_list:
  in_file = open(filename)
  lines = in_file.readlines()
  in_file.close()

  update_file = True
  # look for license separator
  if sep in lines:
    # this file has a license text
    # first, find where it is
    sep_location = lines.index(sep)
    # next, compare if the license text matches the text to add,
    # if so, there is nothing to do. If not, remove the old text
    # to have the new text added
    if license_lines == lines[:sep_location]:
      update_file = False
    else:
      lines = lines[sep_location + 1:]
  
  if (update_file):
    lines = license_lines + [sep,] + lines    
    out_file = open(filename, 'wb')
    for line in lines:
      out_file.write(line)
    out_file.close()

