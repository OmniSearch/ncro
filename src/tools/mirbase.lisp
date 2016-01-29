;; class to keep indexes and entry information.

(defclass mirbase ()
  ((accession2name :accessor accession2name :initform (make-hash-table :test 'equalp))
   (name2accession :accessor name2accession :initform (make-hash-table :test 'equalp))
   (member2family :accessor member2family :initform (make-hash-table :test 'equalp))
   (family2member :accessor family2member :initform (make-hash-table :test 'equalp))
   (entries :accessor entries)
   (problem-families :initform '("MIPF0000783" "MIPF0000773" "MIPF0000153") :allocation :class :accessor problem-families) ; family name matches member name
  ))

;; reads the alias file, which maps miRBase identifiers with labels
(defmethod read-aliases ((m mirbase) &optional (aliases "ncro:src;mirbase;aliases.txt"))
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

(defmethod read-families ((m mirbase) &optional (families "ncro:src;mirbase;miFam.dat"))
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
			  (setf (gethash accession (accession2name m)) name)		    
			  (setf (gethash accession (member2family m)) current-accession))

		     (unless (not (gethash current-id (name2accession m))) (debug "current name already present"))
		     (setf (gethash current-id (name2accession m)) current-accession)

		     (unless (not (gethash current-accession (accession2name m))) (debug "current-accession already registered"))
		     (setf (gethash current-accession (accession2name m)) current-id)))
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

(defmethod read-entries ((m mirbase) &optional (families "ncro:src;mirbase;miRNA.dat"))
  (let ((state :nothing) 
	(current-description nil)
	(current-name nil)
	(current-accession nil)
	(current-longname nil)
	(current-database-references nil))
    (labels ((debug (string)
	       (print string)
	       (print-db state current-name current-accession current-description current-longname current-database-references))
	     (finish-entry ()
	       ;(debug "entry")
	       (setq state :nothing current-members nil current-accession nil current-description nil current-longname nil current-database-references nil))
	     (set-accession (accession)
	       (if current-accession
		   (debug "accession twice in a record")
		   (setq current-accession accession state :accession)
		   ))
	     (set-longname (string)
	       (if current-longname
		   (debug "already a long name")
		   (setq current-longname string state :longname))
	       )
	     (set-description (string)
	       (push string current-description)
	       (setq state :description))
	     (set-database-references (string)
	       (push string current-database-references)
	       (setq state :references))
	     (skip (where line field)
	       (declare (ignore where line field))
	       '(print-db where line field) ))
      (with-open-file (f families :direction :input)
	(loop for line = (read-line f nil :eof)
	   until (eq line :eof)
	   for (field  restline) = (car (all-matches line "(..)\\s*(.*)" 1 2))
	   do
	     (if (eq state :sequence)
		 (if (not (equal field "//"))
		     (skip "sequence" line field)
		     (finish-entry))
		 (cond ((equal field "//")
			(finish-entry))
		       ((member field '("ID" "RA" "RN" "RT" "RL"  "FH" "FT" "RC" "RX" "XX") :test 'equalp) (skip "member" line field))
		       ;;((equal field "ID")
		       ;; (set-id restline))
		       ((equal field "AC")
			(set-accession restline))
		       ((equal field "DE")
			(set-longname restline))
		       ((equal field "DR")
			(set-database-references restline))
		       ((equal field "CC")
			(set-description restline))
		       ((equal field "SQ")
			(setq state :sequence)
			(skip "SQ" line field)
			)
		       (t (debug (format nil "unexpected line: '~a'" line)))))
	   finally (if (and current-accession (not (eq state :nothing)))
		       (finish-entry)))))))



#|
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
