;TODO write "for x over x*" generator to iterate over subsequences.
(defpackage :infr/utils
  (:shadow :mean)
  (:shadowing-import-from :num-utils :sum)
  (:use :cl :alexandria :iterate :num-utils)
  (:export :elmt :mean :vector-sum :vector-mean))

(in-package :infr/utils)

(defun elmt (sequence index)
  "elt with negtive indices. (-1 for last elt. etc.)"
  (elt sequence (mod index (length sequence))))

(defmacro-clause (mean x &optional with-weight (w 1) into var)
  (with-gensyms (num div)
    `(progn (reducing ,x by (lambda (acc x) (+ acc (* ,w x))) initial-value 0 into ,num)
            (summing ,w into ,div)
            ,(if var `(for ,var = (/ ,num ,div))
                     `(finally (return (/ ,num ,div)))))))

(defmacro-clause (vector-sum v &optional into var) `(reducing ,v by #'e+ initial-value 0 into ,var))
(defmacro-clause (vector-multiply v &optional into var) `(reducing ,v by #'e* initial-value 0 into ,var))

(defmacro-clause (vector-mean v &optional with-weight (w 1) into var)
  (with-gensyms (num div)
    `(progn (reducing ,v by (lambda (acc v) (e+ acc (e* ,w v))) initial-value 0 into ,num)
            (summing ,w into ,div)
            ,(if var `(for ,var = (e/ ,num ,div))
                     `(finally (return (e/ ,num ,div)))))))
