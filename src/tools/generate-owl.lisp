(defmethod species-for-prefix ((m mirbase) prefix)
  (second (assoc prefix '(("hsa" !obo:NCBITaxon_9606)
			  ("mmu" !obo:NCBITaxon_10090))
		 :test 'equal)))
  
(defmethod ncro-uri ((m mirbase) entry)
  (make-uri (format nil "http://example.org/~a" (getf entry :accession))))

(defmethod generate-mirbase-owl ((m mirbase))
  (maphash 
   (lambda(accession entry)
     (let ((term-uri (ncro-uri m entry))
	   (only-in-taxon !obo:RO_0002160)
	   (alternative-term !IAO_0000118)
	   (term-species (species-for-prefix m (subseq (car (getf entry :name)) 0 3))))
       (when term-species 
	 (print `((declaration (class ,term-uri ))
		  (annotation !rdfs:label ,(car (getf entry :name)))
		  ,@(loop for label in (cdr (getf entry :name))
		       collect 
			 (list 'annotation alternative-term label))
		  (annotation !oboinowl:hasDbXref ,accession)
		  (annotation ,alternative-term ,(getf entry :longname))
		  (subclass-of ,term-uri !obo:NCRO_0004001)
		  (subclass-of ,term-uri (object-some-values-from ,only-in-taxon ,term-species ))
		  )) (break))))
   (entries m))
  )
