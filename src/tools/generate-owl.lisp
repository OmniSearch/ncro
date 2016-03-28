;; common terms.
(def-uri-alias "definition" !obo:IAO_0000115)
(def-uri-alias "editornote" !obo:IAO_0000116)
(def-uri-alias "termeditor" !obo:IAO_0000117)
(def-uri-alias "definitionsource" !obo:IAO_0000119)
(def-uri-alias "onlyintaxon" !obo:RO_0002160)
(def-uri-alias "alternativeterm" !obo:IAO_0000118)
(def-uri-alias "derivesfrom" !obo:RO_0001000)
(def-uri-alias "maturemirna" !obo:NCRO_0004001)
(def-uri-alias "disjointkey" !obo:ncro-family-disjoint-key)

(defparameter *species-uri-association*
  '(("hsa" !obo:NCBITaxon_9606 "Homo sapiens")
    ("mmu" !obo:NCBITaxon_10090 "Mus musculus")
    ("rno" !obo:NCBITaxon_10116 "Rattus norvegicus")
    ("gga" !obo:NCBITaxon_9031 "Gallus gallus")
    ("cel" !obo:NCBITaxon_6239 "C.elegans")
    ("cfa" !obo:NCBITaxon_9615 "Canis lupus familiaris")
    ("dre" !obo:NCBITaxon_7955 "Danio rerio")
    ("dan" !obo:NCBITaxon_7217 "Drosophila ananassae")
    ("der" !obo:NCBITaxon_7220 "Drosophila erectus")
    ("dgr" !obo:NCBITaxon_7222 "Drosophila grimshawi")
    ("dme" !obo:NCBITaxon_7227 "Drosophila melanogaster")
    ("dmo" !obo:NCBITaxon_7230 "Drosophila mojavensis")
    ("dpe" !obo:NCBITaxon_7234 "Drosophila persimilis") 
    ("dps" !obo:NCBITaxon_7237 "Drosophila pseudoobscura")
    ("dse" !obo:NCBITaxon_7238 "Drosophila sechellia")
    ("dsi" !obo:NCBITaxon_7240 "Drosophila simulans")
    ("dvi" !obo:NCBITaxon_7244 "Drosophila virilis")
    ("dwi" !obo:NCBITaxon_7260 "Drosophila willistoni")
    ("dya" !obo:NCBITaxon_7245 "Drosophila yakuba")))

(defmethod species-for-prefix ((m mirbase) prefix)
  (second (assoc prefix *species-uri-association*
		 :test 'equalp)))

(defmethod species-name-for-prefix ((m mirbase) prefix)
  (third (assoc prefix *species-uri-association*
		 :test 'equalp)))
  
(defmethod uri-for-accession ((m mirbase) accession)
  (make-uri (format nil "http://purl.obolibrary.org/NCRO_~a" accession)))

#|
The issue here is how to align what is wanted (as best I can guess)
from NCRO and the work that has been done already, primary the
associate with miRBase. Currently the dxref is to miRBase records that
appear to denote pre-miRNA. I conclude this because the sequences
given include the intact stem-loop and modified caps but not the
poly-a and because miRBase has different entries for mature miRNA (the
MIMATxxx records).

NCRO initially mixes these records. We have, for example,

hsa-miR-125b-5p (mature - MIMAT) MIMAT0000423 from hsa-miR-125b-1 (this should really be called  hsa-miR-125b-1-5p )
hsa-miR-125b-1-3p (mature - MIMAT) MIMAT0004592 from hsa-miR-125b-1
hsa-miR-125b-2-3p (mature - MIMAT) MIMAT0004603 from hsa-miR-125b-2
hsa-miR-125b-1 (stem loop - MI0) MI0000446
hsa-miR-125b-2 (stem loop) - MI0) MI0000470


It also has
mir-10 MIPF0000033 - the family that includes the above

mir-10 is related to , hsa-miR-125b-1 and hsa-miR-125b-2 via
is_about_grouped_miRNA

The mature forms are not related to the stem loops in the original.

Since both stem loops and mature are present in NCRO (albeit not well
distinguished) we should include both, with a derived-from link from
mature to pre. (derivation is transitive)

Now families. The question is whether the families are computational
artifacts or evolutionary relationship. From reading the paper it
seems that they are intended to be evolutionary, even though they may
be determined to be such by sequence alignment. So we will have
families be superclasses (ala cross-species in PRO) with evidence
coming from the algorithm application.

Finally, there is the question of whether to represent the sequences.
Since the mature miRNA's are unambiguous it seems that might be a good
idea. What is the relation to the molecules? Concretization? The
sequences are individuals.

What should the definitions be?

The gene products should mention the gene. We don't have an ontology
of genes. We'll name the gene symbol when we have one.  The final
products are the mature sequences. We can name them in the
definition. We do have the chromosome positions - we could use that.

The stem loop definitions can either include or not include the sequence.q

The mature miRNA should be identified by the sequence which we can
name in the definition.

The families can be named by the human gene - blah or orthologous to
blah.
|#

;; human, mouse, rat, c-elegams, dros, dog, cow, chicken, zfin
(defmethod human-or-model-in-human-family ((m mirbase) entry)
  (let* ((name (getf entry :name)))
    (or (and (#"matches" name "hsa-.*") entry)
	(and (#"matches" name "^(mmu|rno|d|cel|cfa|gga|dre)-.*")
	     (let ((fam (gethash (#"replaceAll" (getf entry :accession) ";" "") (member2family m))))
	       (find-if (lambda(el) (#"matches" (getf el :name) "hsa.*")) 
			(mapcar (lambda (el) (gethash (concatenate 'string el ";") (entries m)))
				(gethash fam (family2member m)))))))))

(defmethod human-entry ((m mirbase) entry)
  (let* ((name (getf entry :name)))
    (#"matches" name "hsa-.*")))

;; print out and count the mirna to include
(defmethod debug-which-to-include ((m mirbase))
  (let ((count 0))(maphash 
   (lambda(accession entry) ;; iterate over stem loops
     (when (human-or-model-in-human-family m entry)
       (print (getf entry :name)) (incf count)))
     (entries m)) count))

;; FIXME This could be called more than once because of more than one human entry for a family
(defmethod axioms-for-family ((m mirbase) human-entry family-accession)
  (let ((subject (uri-for-accession m family-accession))
	(mirna-gene-product !obo:NCRO_0004019)
	(human-members (map 'list 'car
			  (remove-if-not 
			   (lambda(el)(#"matches" (car el) "hsa.*"))
			   (loop for member in (gethash family-accession (family2member m) )
			      collect (gethash member (accession2name m)))))))
    (list*
     `(declaration (class ,subject))
     `(subclassof ,subject ,mirna-gene-product)
     `(annotation-assertion !foaf:page ,subject ,(make-uri (format nil "http://www.mirbase.org/cgi-bin/mirna_summary.pl?fam=~a" family-accession)))
     `(annotation-assertion !rdfs:label ,subject ,(format nil "miRNA family ~a" (gethash family-accession (accession2name m))))

     `(annotation-assertion
       !definition ,subject
       ,(format nil "An evolutionary conserved group of miRNA gene products including the human ~{~a~^, ~}, their precursors and their orthologs and paralogs" human-members))
     `(annotation-assertion !definitionsource ,subject "miRBase")
     `(annotation-assertion !alternativeterm ,subject ,family-accession)
     `(subclassof ,subject (data-has-value !disjointkey ,family-accession))
     (loop for member in (gethash family-accession (family2member m))
	for entry = (gethash (concatenate 'string member ";") (entries m))
	when (human-or-model-in-human-family m entry)
	collect
	  `(subclassof ,(uri-for-accession m member) ,subject)))))

(defmethod axioms-for-pre-mirna ((m mirbase) entry human-entry)
  (let ((subject (uri-for-accession m (#"replaceAll" (getf entry :accession) ";" "")))
	(gene (gethash (#"replaceAll" (getf human-entry :accession) ";" "") (accession2gene m)))
	(pre-mirna !obo:NCRO_0004020))
    (list
     `(declaration (class ,subject))
     `(annotation-assertion !rdfs:label ,subject ,(getf entry :longname))
     (if (getf entry :description)
	 `(annotation-assertion !editornote ,subject ,(getf entry :description)))
     `(subclassof ,subject ,pre-mirna)
     `(annotation-assertion !foaf:page ,subject ,(make-uri (format nil "http://www.mirbase.org/cgi-bin/mirna_entry.pl?acc=~a" (getf entry :accession))))
     `(annotation-assertion !definition ,subject
			    ,(format nil "A pre-miRNA processed from the primary transcript~a of human gene ~a" (if (eq entry human-entry) "" (format nil " of a ~a ortholog" (species-for-prefix m (subseq (getf entry :name) 0 3)))) (first gene)))
     `(subclass-of ,subject (object-some-values-from !onlyintaxon ,(species-for-prefix m (subseq (getf entry :name) 0 3))))
     `(annotation-assertion !alternativeterm ,subject ,(#"replaceAll" (getf entry :accession) ";" ""))
     )))


(defmethod axioms-for-matures ((m mirbase) entry)
  (let ((pre-mirna-uri  (uri-for-accession m (#"replaceAll" (getf entry :accession) ";" "")))
	(family (gethash (#"replaceAll" (getf entry :accession) ";" "") (member2family m))))
    (loop for mature in (getf entry :matures)
       for sequence = (second (assoc :sequence mature))
       for name = (second (assoc :product mature))
       for accession = (second (assoc :accession mature))
       for subject = (uri-for-accession m accession)
       append
	 (list
	  `(declaration (class ,subject))
	  `(annotation-assertion !rdfs:label ,subject  ,name)
	  `(subclassof ,subject !maturemirna)
	  (when family `(subclassof ,subject ,(uri-for-accession m family)))
	  `(subclassof ,subject (object-some-values-from !derivesfrom ,pre-mirna-uri))
	  `(annotation-assertion !foaf:page ,subject ,(make-uri (format nil "http://www.mirbase.org/cgi-bin/mature.pl?mature_acc=~a" accession)))
	  `(annotation-assertion !definition ,subject
				 ,(format nil "A mature miRNA expressed in ~a, with sequence ~a" (species-for-prefix m (subseq (getf entry :name) 0 3)) sequence))
	  `(subclass-of ,subject (object-some-values-from !onlyintaxon ,(species-for-prefix m (subseq (getf entry :name) 0 3))))
	  `(annotation-assertion !alternativeterm ,subject ,accession)))))


(defmethod generate-mirbase-owl ((m mirbase))
  (with-ontology mirna (:collecting t :ontology-iri "http://purl.obolibrary.org/obo/ncro/dev/ncro-mirna.owl") ;; :only-return-axioms t :also-return-axioms t)
      ((let ((seen (make-hash-table :test 'equalp)))
	 (asq (imports !obo:ncro.owl)
	      (declaration (data-property !disjointkey))
	      (functional-data-property !disjointkey)
	      (annotation !definition !disjointkey "A computational artifact that allows a more compact way of ensuring that the miRNA families are considered disjoint"))
	 (maphash 
	  (lambda(accession entry) ;; iterate over stem loops
	    (let ((human-entry (human-or-model-in-human-family m entry)))
	      (when human-entry
		(as (axioms-for-pre-mirna m entry human-entry))
		(as (axioms-for-matures m entry))
		(when (human-entry m entry)
		  (let ((fam (gethash (#"replaceAll" (getf entry :accession) ";" "") (member2family m))))
		    (when fam
		      (unless (gethash fam seen)
			(as (axioms-for-family m entry (gethash (#"replaceAll" (getf entry :accession) ";" "") (member2family m))))
			(setf (gethash fam seen) t)
			)))))))
	  (entries m))))
    (write-rdfxml mirna "~/Desktop/ncro-mirna.owl")
    ))

(defun species-presentation-name (entry )

    (species-for-prefix m (subseq (getf entry :name) 0 3)))
