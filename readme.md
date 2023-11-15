# IMGT Germlines

Two projects:
* generate 
* library

Generate generates the binary germlines and rust code needed to interface with them. The library can then be pushed to crates.io while not overflowing the storage limit. If you want to generate newer germlines from an update database you need to put the `imgt.dat.Z` file in the `data` directory (this can be downloaded from `https://www.imgt.org/download/LIGM-DB/imgt.dat.Z`).

Folders:
* data - put the IMGT data here
* generate - the rust project to generate the germline binary files
* germlines - the generated binary files + rust code to use them
* library - the library which ends up on crates.io
* shared - shared code between the library and generate binary