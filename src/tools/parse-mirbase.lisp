;; class to keep indexes and entry information.

(defclass mirbase ()
  ((accession2name :accessor accession2name :initform (make-hash-table :test 'equalp))
   (name2accession :accessor name2accession :initform (make-hash-table :test 'equalp))
   (member2family :accessor member2family :initform (make-hash-table :test 'equalp))
   (family2member :accessor family2member :initform (make-hash-table :test 'equalp))
   (accession2gene :accessor accession2gene :initform (make-hash-table :test 'equalp))
   (accession2sequence :accessor accession2sequence :initform (make-hash-table :test 'equalp))
   (entries :accessor entries :initform (make-hash-table :test 'equalp))
   (problem-families :initform '("MIPF0000783" "MIPF0000773" "MIPF0000153") :allocation :class :accessor problem-families) ; family name matches member name
   (ncro :allocation :class)
  ))

;; reads the alias file, which maps miRBase identifiers with labels
(defmethod read-aliases ((m mirbase) &optional (aliases "ncro:src;mirbase;21;aliases.txt"))
  (with-open-file (f aliases :direction :input)
    (loop with 2name = (accession2name m)
       with 2accession = (name2accession m)
       for line = (read-line f nil :eof)
       until (eq line :eof)
       for (mirbaseid names) = (split-at-char line #\tab)
       do (setq names (split-at-char names #\;))
	 (setf (gethash mirbaseid 2name) names)
	 (loop for name in names do
	      (pushnew mirbaseid (gethash name 2accession) :test 'equalp))))
  m)

;; check that a list of requests (one name per request) are actually names in mirbase
(defmethod check-request ((m mirbase) &optional (requests "ncro:src;requests.txt"))
  (let ((requests 
	 (with-open-file (f requests :direction :input)
	   (loop for line = (read-line f nil :eof)
	      until (eq line :eof)
	      collect line))))
    (loop for request in requests
       when (not (gethash request (name2accession m)))
       do (format t "~a not in mirbase~%" request))))

;; read the families file and build tables mapping members to families and vice versa

(defmethod read-families ((m mirbase) &optional (families "ncro:src;mirbase;21;miFam.dat"))
  (let ((state :nothing) (current-members nil) (current-id nil) (current-accession nil))
    (labels ((debug (string)
	       (print string)
	       (print-db state current-members current-id current-accession))
	     (finish-entry ()
	       (if (member current-accession (problem-families m) :test 'equalp)
		   (format t "skipping family ~a ~a~%" current-accession current-id)
		   (progn
		     (unless (and current-members current-id current-accession) (debug "something missing!"))
		     (setf (gethash current-accession (family2member m)) (mapcar 'first current-members))
		     (loop for (accession . name) in current-members do
			  (unless (not (gethash accession (member2family m))) (debug "Member already has a family"))
			  (setf (gethash name (name2accession m)) accession)
			  (pushnew name (gethash accession (accession2name m)) :test 'equal)
			  (setf (gethash accession (member2family m)) current-accession))

		     (unless (not (gethash current-id (name2accession m))) (debug "current name already present"))
		     (setf (gethash current-id (name2accession m)) current-accession)

		     (unless (not (gethash current-accession (accession2name m))) (debug "current-accession already registered"))
		     (pushnew current-id (gethash current-accession (accession2name m)) :test 'equal)))
	       (setq state :nothing current-members nil current-id nil current-accession nil)
	       )
	     (set-id (id)
	       (assert (null current-id) (state current-members current-id current-accession)
		       "ID field appears twice in a record")
	       (setq current-id id state :gather)
	       )
	     (set-accession (accession)
	       (assert (null current-accession) (state current-members current-id current-accession)
		       "accession field appears twice in a record")
	       (setq current-accession accession state :gather)
	       )
	     (add-member (accession name)
	       (assert (not (eq state :nothing)) (state) "should have had an id by now")
	       (push (cons accession name) current-members)
	       (setq state :gather)
	       ))
      (with-open-file (f families :direction :input)
	(loop for line = (read-line f nil :eof)
	   with state = :nothing
	   with current = nil
	   until (eq line :eof)
	   for (field . rest) = (split-at-regex line "\\s+")
	   do
	     (cond ((equal field "//")
		     (finish-entry))
		    ((equal field "ID")
		     (set-id (first rest)))
		    ((equal field "AC")
		     (set-accession (first rest)))
		    ((equal field "MI")
		     (add-member (first rest) (second rest)))
		    (t (debug (format nil "unexpected line: '~a'" line))))
	   finally (if (and (equal current nil) (not (eq state :nothing))) (finish-entry)))))))

;; read the entries in miRNA.dat, only paying attention to accession, long name, description, and database references

(defmethod create-entry ((m mirbase) accession name longname description sequence matures hgnc)
  (setf (gethash accession (entries m))
	`(:accession ,accession :gene (:hgnc ,@hgnc) :name ,name :longname ,longname :description ,(#"replaceAll" (format nil "~{~a~^ ~}" (reverse description)) "\\s+" " ") :sequence ,sequence :matures ,matures))
  )
			 
(defmethod read-entries ((m mirbase) &optional (families "ncro:src;mirbase;21;miRNA.dat"))
  (let ((state :nothing) 
	(current-description nil) (current-name nil) (current-accession nil)
	(current-longname nil) (current-database-references nil) (current-sequence nil)	(current-matures nil)
	(current-hgnc nil))
    (labels ((debug (string)
	       (print string)
	       (print-db state current-name current-accession current-description current-longname current-database-references))
	     (finish-entry ()
	       (create-entry m current-accession current-name current-longname
			     current-description current-sequence current-matures current-hgnc)
	       (setq state :nothing current-members nil current-sequence nil
		     current-name nil current-accession nil current-description nil
		     current-longname nil current-database-references nil current-matures nil
		     current-hgnc nil)
	       )
	     (skip (where line field)
	       (declare (ignore where line field))
	       '(print-db where line field) )
	     )
      (with-open-file (f families :direction :input)
	(loop with peeked
	   for line = (or peeked (read-line f nil :eof))
	   until (eq line :eof)
	   for (field  restline) = (car (all-matches line "(..)\\s*(.*)" 1 2))
	   do
	     (setq peeked nil)
	     (if (eq state :sequence)
		 (if (not (equal field "//"))
		     (skip "sequence" line field)
		     (finish-entry))
		 (cond ((equal field "//")
			(finish-entry))
		       ((member field '("RA" "RN" "RT" "RL"  "FH" "RC" "RX" "XX") :test 'equalp) (skip "member" line field))
		       ((equal field "ID")
			(if current-name
			    (debug "already a current name")
			    (let ((parts (car (all-matches restline "(.*?)\\s+(.*?); (.*?); (.*?); (.*?)." 1))))
			      (setq current-name (car parts) state :name))))
		       ((equal field "AC")
			(if current-accession
			    (debug "accession twice in a record")
			    (setq current-accession restline state :accession)
			    ))
		       ((equal field "DE")
			(if current-longname
			    (debug "already a long name")
			    (setq current-longname restline state :longname)))
		       ((equal field "DR")
			;; DR   HGNC; 35259; MIR1178.
			(let ((hgnc (car (all-matches restline "(HGNC|ENTREZGENE);\\s*(\\d+);\\s*(.*)\\." 2 3))));; FIXME - not all gene names are from HGNC
			  (and hgnc (setq current-hgnc hgnc))
			  (push restline current-database-references)))
		       ((equal field "CC")
			(push restline current-description))
		       ((equal field "SQ")
			(setq state :sequence)
			(loop for line = (read-line f) until (equal line "//") collect line into lines
			   finally
			     (let ((concat (apply 'concatenate 'string lines)))
			       (setq concat (#"replaceAll" concat "\\d|\\s" ""))
			       (unless (#"matches" concat "^[ucgarywsmkhbvdn]+$") (format t "misread sequence for ~a: ~a~%" current-accession concat))
			       (setq current-sequence concat)
			       (finish-entry))))
			 
		       ((equal field "FT")
			(if (#"matches" line "^FT\\s+(modified_base|Key).*")
			    (read-line f)
			    (destructuring-bind (&optional from to)
				(car (all-matches restline "miRNA\\s+(\\d+)\\.\\.(\\d+)" 1 2))
			      (unless (and from to) (break))
			      (loop for line = (read-line f nil :eof)
				 with props
				 until (or (#"matches" line "^FT\\s+miRNA.*") (not (#"matches" line "^FT.*")))
				 do (destructuring-bind (&optional key value)
					(car (all-matches line "FT\\s+/(\\w+)=\"([^\"]+)" 1 2))
				      (and key
					   (push (list (intern (string-upcase (string key)) 'keyword) value) props))
				      )
				 finally
				   (progn
				     (push `((:from ,from) (:to ,to) ,@props) current-matures)
				     (setq peeked line))))))
		       (t (debug (format nil "unexpected line: '~a'" line)))))
	   finally (if (and current-accession (not (eq state :nothing)))
		       (finish-entry))
	     )))))

(defmethod check-ncro-mirna-against-mirbase ((m mirbase) &optional (ontology-file "ncro:src;ontology;ncro.owl"))
  (unless (and (slot-boundp m 'ncro) (slot-value m 'ncro))
    (setf (slot-value m 'ncro)
	  (load-ontology ontology-file)))
  (let ((ontology (slot-value m 'ncro))
	(mirna-class !obo:NCRO_0004001))
    (let ((classes (descendants  mirna-class ontology)))
      (loop for class in classes
	 for label = (car (rdfs-label class ontology))
	 for dbxref = (second (first (entity-annotations class ontology !oboinowl:hasDbXref)))
	 for mirbase-accession = (caar (all-matches dbxref "miRBase:(.*)" 1))
	 do
	   (when dbxref
	     (let ((lookup-names (gethash mirbase-accession (accession2name m))))
	       (when (not (member label lookup-names :test 'equal))
		 (format t "~a	~{~a~^;~}		~a	~a~%" mirbase-accession lookup-names label (gethash label (name2accession m))))))))))

;; chr1	.	miRNA_primary_transcript	17369	17436	.	-	.	ID=MI0022705;Alias=MI0022705;Name=hsa-mir-6859-1

(defmethod read-human-gene-positions ((m mirbase) &optional (file "ncro:src;mirbase;21;genomes;hsa.gff3"))
  (with-open-file (f file :direction :input)
    ;; skip comments
    (loop for peek = (peek-char nil f)
       until (not (char= peek #\#))
       do (read-line f))
    (loop for line = (read-line f nil :eof)
       until (eq line :eof)
       for (chromosome nil nil from to nil nil nil which) = (split-at-char line #\tab)
       for id = (caar (all-matches which "ID=(MI\\d+);" 1))
       do (setf (gethash id (accession2gene m)) (list (third (getf (gethash (concatenate 'string id ";") (entries m)) :gene)) chromosome from to)))))

;; Using the FT annotations and the stem loop sequence, compute the mature sequence
;; So far spot-checked - should validate with independent source
(defmethod compute-mature-sequence ((m mirbase))
  (maphash (lambda(accession entry)
	     (loop for mature in (getf entry :matures)
		for from = (second (assoc :from mature))
		for to = (second (assoc :to mature))
		do
		  (cond ((and from (not to)) (and to (not from))
			 (warn "~a mature sequence not delimited well; ~a" accession mature))
			((not (or from to)) (print-db entry mature from to) (break))
			(t (nconc mature `((:sequence ,(subseq  (getf entry :sequence) (1- (read-from-string from)) (read-from-string to)))))))))
	   (entries m)))
		  
;; check if any of the sequence letters are other than the standard
;; uagc (shouldn't be, even though the stem loop sequence may be
;; ambiguous

(defmethod ambiguous-mature-sequences ((m mirbase))
  (maphash (lambda(accession entry)
	     (loop for mature in (getf entry :matures)
		for sequence = (second (assoc :sequence mature))
		when (not (#"matches" sequence "^[uagc]+$"))
		  do (warn "ambiguous mature sequence for ~a: ~a" accession sequence)))
	   (entries m)))

(defun parse-mirbase ()
  (let ((m (make-instance 'mirbase)))
    (read-aliases m)
    (read-entries m)
    (read-families m)
    (compute-mature-sequence m)
    (read-human-gene-positions m)
    m))

#|

CiteXplore was withdrawn from service on 15 February 2013, replaced by Europe PubMed Central
https://github.com/ebi-chebi/ChEBI/issues/3134

Family level
 Gene level (species specific - gene product of .... )
   sequence level (with stem loop)
   sequence level (after stem loop removed) (derives from above)
   mature (single stranded 20-24nt) (derives from above)
    5' 
    3' 

When two mature microRNAs originate from opposite arms of the same
pre-miRNA and are found in roughly similar amounts, they are denoted
with a -3p or -5p suffix. (In the past, this distinction was also made
with 's' (sense) and 'as' (antisense)). However, the mature microRNA
found from one arm of the hairpin is usually much more abundant than
that found from the other arm,[2] in which case, an asterisk following
the name indicates the mature species found at low levels from the
opposite arm of a hairpin. For example, miR-124 and miR-124* share a
pre-miRNA hairpin, but much more miR-124 is found in the cell

http://www.chem.qmul.ac.uk/iupac/bibliog/white.html
http://www.chem.qmul.ac.uk/iubmb/misc/naseq.html
3. Allocation of symbols
3.1. Guanine. adenine, thymine, cytosine: G,A,T,C
3.2. Purine (adenine or guanine): R
3.3. Pyrimidine (thymine or cytosine): Y
3.4. Adenine or thymine: W
3.5. Guanine or cytosine: S
3.6. Adenine or cytosine: M
3.7. Guanine or thymine: K
3.8. Adenine or thymine or cytosine: H
3.9. Guanine or cytosine or thymine: B
3.10. Guanine or adenine or cytosine: V
3.11. Guanine or adenine or thymine: D
3.12. Guanine or adenine or thymine or cytosine: N

FT   miRNA           17..38
FT                   /accession="MIMAT0000001"
FT                   /product="cel-let-7-5p"
FT                   /evidence=experimental
FT                   /experiment="cloned [1-3], Northern [1], PCR [4], 454 [5],
FT                   Illumina [6], CLIPseq [7]"


miFAM.dat - entries
AC - Accession 
ID - ID
MI - mirna name
//

			  

//

SQ
XX

miRNA.dat.zip - entries
miFam.dat species-specific to family

XX - blank line
ID - species-specific name
AC - mirbase accession
DE - long name
RN - reference number
RX - reference identifier
RA - authors
RT - title
RL - journal
DR - database reference
CC - description
FH - feature header
FT - feature body
RC - some other sort of article reference
SQ - sequence. continuation lines don't have two-letter field name 

Sample record

ID   cel-let-7         standard; RNA; CEL; 99 BP.
XX
AC   MI0000001;
XX
DE   Caenorhabditis elegans let-7 stem-loop
XX
RN   [1]
RX   PUBMED; 11679671.
RA   Lau NC, Lim LP, Weinstein EG, Bartel DP;
RT   "An abundant class of tiny RNAs with probable regulatory roles in
RT   Caenorhabditis elegans";
RL   Science. 294:858-862(2001).
XX
RN   [2]
RX   PUBMED; 12672692.
RA   Lim LP, Lau NC, Weinstein EG, Abdelhakim A, Yekta S, Rhoades MW, Burge CB,
RA   Bartel DP;
RT   "The microRNAs of Caenorhabditis elegans";
RL   Genes Dev. 17:991-1008(2003).
XX
RN   [3]
RX   PUBMED; 12747828.
RA   Ambros V, Lee RC, Lavanway A, Williams PT, Jewell D;
RT   "MicroRNAs and other tiny endogenous RNAs in C. elegans";
RL   Curr Biol. 13:807-818(2003).
XX
RN   [4]
RX   PUBMED; 12769849.
RA   Grad Y, Aach J, Hayes GD, Reinhart BJ, Church GM, Ruvkun G, Kim J;
RT   "Computational and experimental identification of C. elegans microRNAs";
RL   Mol Cell. 11:1253-1263(2003).
XX
RN   [5]
RX   PUBMED; 17174894.
RA   Ruby JG, Jan C, Player C, Axtell MJ, Lee W, Nusbaum C, Ge H, Bartel DP;
RT   "Large-scale sequencing reveals 21U-RNAs and additional microRNAs and
RT   endogenous siRNAs in C. elegans";
RL   Cell. 127:1193-1207(2006).
XX
RN   [6]
RX   PUBMED; 19460142.
RA   Kato M, de Lencastre A, Pincus Z, Slack FJ;
RT   "Dynamic expression of small non-coding RNAs, including novel microRNAs
RT   and piRNAs/21U-RNAs, during Caenorhabditis elegans development";
RL   Genome Biol. 10:R54(2009).
XX
RN   [7]
RX   PUBMED; 20062054.
RA   Zisoulis DG, Lovci MT, Wilbert ML, Hutt KR, Liang TY, Pasquinelli AE, Yeo
RA   GW;
RT   "Comprehensive discovery of endogenous Argonaute binding sites in
RT   Caenorhabditis elegans";
RL   Nat Struct Mol Biol. 17:173-179(2010).
XX
DR   RFAM; RF00027; let-7.
DR   WORMBASE; C05G5/12462-12364; .
XX
|#
