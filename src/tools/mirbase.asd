(in-package :asdf)

(setf (logical-pathname-translations "ncro")
      `(("**;*.*" ,(make-pathname :directory (append (butlast (pathname-directory *load-pathname*) 2)
						     '(:wild-inferiors))
				  :name :wild
				  :type :wild
				  :defaults *load-pathname*))))

(defsystem :ncro
  :name "NCRO Tools"
  :author "Alan Ruttenberg"
  :license "BSD"
  :components
  ((:file "parse-mirbase")
   (:file "generate-owl")))

