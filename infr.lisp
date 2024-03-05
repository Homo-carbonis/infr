(require 'lisp-stat)
(require 'utils)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr
  (:shadowing-import-from :lisp-stat :product :next :sum :generate :mean)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :iter-utils :lisp-stat :iter)
  (:export :generate-markov-chain :estimate-parameters :likelihood :marginal :posterior :posteriors))

(in-package :infr)

(defun generate-markov-chain (f n &key (xi #(0d0)))
  (iter (with offset = (length xi))
        (with x = (make-array offset
                              :initial-contents xi
                              :element-type 'double-float
                              :fill-pointer offset))
        (initially (adjust-array x n))
        (for i from offset below n)
        (vector-push (draw (funcall f x)) x)
        (finally (return x))))

; Possibly we should weight the mean by magnitude of y in some way?
(defun estimate-parameters (model prior y &key (offset 2) (sample-count 100))
  "Find a Monte Carlo estimate of the parameter of a model.  We draw a large
   number of samples from 'prior. Then for each sampled value of beta we find
   P(y|model,beta). Finally we take the mean of the samples using P(y|model,beta)
   as weights, to find an expected value for beta."
  (iter (with beta* = (generate (lambda () (each #'draw prior)) sample-count))
        (for i from offset below (length y)) 
        (setf (fill-pointer y) i) ; we pass y up to but not including element i to the model function.
        (vector-mean
          (iter (for beta in-vector beta*)
                (vector-mean beta weight (pdf (funcall model beta y) (aref y i)))))
        (finally (incf (fill-pointer y))))); Restore fill pointer before returning

(defun log-likelihood (model beta y &key (offset 2))
  "Find log P(y|model,beta)."
  (iter (for i from offset below (length y))
        (setf (fill-pointer y) i)
        (summing (log-pdf (funcall model beta y) (aref y i)))
        (finally (incf (fill-pointer y)))))

(defun log-marginal (model prior y &key (offset 2) (sample-count 100))
  "Find log P(y|model)."
  (let* ((n (length prior))
         (beta* (make-array `(,sample-count ,n)
                            :element-type 'double-float))
         (maxima (make-array n :element-type 'double-float))
         (minima (make-array n :element-type 'double-float))
         (volume 0d0))
    (iter (for i below n)
          (iter (for j below sample-count)
                (for beta = (draw (aref prior i)))
                (maxf (aref maxima i) beta)
                (minf (aref minima i) beta)
                (setf (aref beta* j i) beta)))
    (setf volume (product (e- maxima minima)))
    (iter (for beta in-vector (split beta* 1))
          (for p = (pdf* prior beta))
          (for likelihood = (log-likelihood model beta y :offset offset))
          (summing (* p (+ (log p) likelihood)) into num)
          (summing p into div)
          (finally (return (+ (log volume) (/ num div)))))))

(defun pdf* (rv* x*)
  (iter (for rv in-vector rv*)
        (for x in-vector x*)
        (multiplying (pdf rv x))))

(defun log-pdf* (rv* x*)
  (iter (for rv in-vector rv*)
        (for x in-vector x*)
        (summing (log-pdf rv x))))

(defun log-posterior (model prior beta y &key (offset 2) (sample-count 100))
  (let ((likelihood (log-likelihood model beta y))
        (prior (log-pdf* prior beta))
        (marginal (log-marginal model prior y)))
    (print likelihood)
    (print prior)
    (print marginal)
    (exp (+ likelihood prior 
            (- marginal)))))

(defun posteriors (models priors y &key (offset 2))
  "Assign relative probabilities to each of models"
  (let* ((joint* (each (lambda (m p) (* p (likelihood m y :offset offset))) models priors))
         (marginal (sum joint*)))
    (e/ joint* marginal)))
