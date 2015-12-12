#
# This file is part of the DMS molecular simulation package.
#
# Copyright (c) 2014,2015 by the DMS development team.
#
# DMS is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1
# of the License, or (at your option) any later version.
#
# DMS is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with DMS; if not, see
# http://www.gnu.org/licenses, or write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA.
#
# If you want to redistribute modifications to DMS, please
# consider that scientific software is very special. Version
# control is crucial - bugs must be traceable. We will be happy to
# consider code for inclusion in the official distribution, but
# derived work must not be called official DMS. Details are found
# in the README & COPYING files.
#
# To help us fund DMS development, we humbly ask that you cite
# the research papers on the package.

#TODO: add check that source doesn't contain any untracked files
if(CPACK_SOURCE_PACKAGE_FILE_NAME) #building source package
    get_filename_component(CMAKE_BINARY_DIR ${CPACK_OUTPUT_CONFIG_FILE} PATH)
    if (NOT EXISTS "${CMAKE_BINARY_DIR}/share/man/man1/gmx-view.1" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/INSTALL" OR
        NOT EXISTS "${CMAKE_BINARY_DIR}/share/html/final/online.html")
        message(FATAL_ERROR
            "To create a complete source package all man and HTML pages need "
            "to be generated, and the INSTALL file generated. "
            "Run 'make man html' to build these parts. You can also set "
            "GMX_BUILD_HELP=ON to automatically build the HTML parts. "
            "The latter also requires you to execute 'make install-guide'.")
    endif()
endif()
