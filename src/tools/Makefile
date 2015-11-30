# Makefile 
#
# It is intended to be run on a Unix system (Linux, Mac OS X)
#
# WARNING: Makefiles contain significant tab characters!
# Ensure that your editor shows tab characters before editing this file.


### Configuration

imports:../ontology/chebi-imports.owl
	curl -s -F file=@../ontology/ontofox-inputs/chebi.txt -o ../ontology/chebi-imports.owl  http://ontofox.hegroup.org/service.php
