Build artifacts, such as releases, should be created in this directory, *in the release branch*, not in the master branch.

The process is

- Create a branch for the new release. Versions are named by date, so the branch would be called release-YYYY-MM-DD

- switch to that branch

- Change all the -imports.owl so that they have ontology URI http://purl.obolibrary.org/obo/ncro/xxx-import.owl and version IRI and version IRI http://purl.obolibrary.org/obo/ncro/YYYY-MM-DD/xxx-import.owl

- Change all the /dev/ URIs in ncro.owl to import by versionIRI instead of /dev

- Before saving, check that "Use XML Entities" is *unchecked* in the "Save" tab of preferences. Create a merged file. I used protege's merge function. Save the merged file in build/ncro-all-in-one.owl

- Fix ontology properties in the merged file. Since the merge lumps all imported ontology properties into the single ontology, I remove those and replace with only the ontology properties in src/ontology/ncro.owl. 

- To the merged file, add a comment like this one:
><rdfs:comment xml:lang="en">This is the release of December 12, 2015. It will always be available at http://purl.obolibrary.org/ncro/2015-12-10/ncro.owl. While current, this ontology is also available at http://purl.obolibrary.org/obo/ncro.owl
> This file is a merge of NCRO and any imported ontologies or terms. You can access the ontology that existed before the merge, with imports intact, at http://purl.obolibrary.org/obo/ncro/YYYY-MM-DD/ncro.owl or, for the current version, http://purl.obolibrary.org/obo/ncro/prebuild/ncro.owl</rdfs:comment>

- Since we have erased evidence of imports, I add a seeAlso for each ontology we've imported from. In the version of 2015-12-10 that looks like:


>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/bfo.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/chebi.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/clo.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/doid.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/go.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/iao.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/obi.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/pr.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/ro.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/so.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/uberon.owl"/>
>        <rdfs:seeAlso rdf:resource="http://purl.obolibrary.org/obo/iao/ontology-metadata.owl"/>


- commit,  push

Next the PURLs need to be updated. You can edit the PURLs in the github interface. Our configuration is [here](https://github.com/OBOFoundry/purl.obolibrary.org/blob/master/config/ncro.yml). 
The config file has comments, but to reiterate: When there is a new release you need to add a prefix/replacement pairs around line 28. You need to update the "products:" entry to point to the raw merged file. Finally you need to adjust the entry "exact: /prebuild/ncro.owl" so that the replacement points to the raw ncro.owl file.

Once you've finished editing, commit the changes. You will be offered to create a new branch. Do so and submit the pull request. The pull request will be handled by someone among the purl.obolibrary.org developers. Once the pull request is merged between 15 and 30 minutes later the purls should be active. Test them. In particular test http://purl.obolibrary.org/obo/ncro.owl, http://purl.obolibrary.org/obo/ncro/prebuild/ncro.owl and http://purl.obolibrary.org/obo/ncro/YYYY-MM-DD/ncro.owl. To test them, load them by URL into protege and making sure there are no exceptions while loading and that content is what you expect.

Alan Ruttenberg
