(require 'lisp-stat)
(require 'utils)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr
  (:shadowing-import-from :lisp-stat :product :next :sum :generate :mean)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :iter-utils :lisp-stat :iter)
  (:export :generate-wiener :generate-chain :generate-brownian-chain :estimate-parameters :likelihood :posteriors))

(in-package :infr)

(defun generate-chain (f n &key (xi #(0d0)))
  (iter (with offset = (length xi))
        (with x = (make-array offset
                              :initial-contents xi
                              :element-type 'double-float
                              :fill-pointer offset))
        (initially (adjust-array x n))
        (for i from offset below n)
        (vector-push (funcall f x) x)
        (finally (return x))))

(defun generate-brownian-chain (f n &key (xi #(0d0)) (var 1d0))
  (generate-chain (lambda (x)
                     (+ (funcall f x) (draw (r-normal 0d0 var))))
                  n :xi xi))


; Possibly we should weight the mean by magnitude of y in some way?
(defun estimate-parameters (model prior y &key (offset 1) (beta-sample-count 100) (y-sample-count 100) )
  "Find a Monte Carlo estimate of the parameter of a model.
   We draw a large number of samples from 'prior. Then for each sampled value
   of beta and each value of y generate samples from the model and compute the
   average squared difference from y. These are then used as weights to
   calculate an expected value for the parameter beta."
  (iter (with beta* = (generate (lambda () (each #'draw prior)) beta-sample-count))
        (for i from offset below (length y)) 
        (setf (fill-pointer y) i) ; we pass y up to but not including element i to the model function.
        (vector-mean
          (iter (for beta in-vector beta*)
                (for dev = (iter (for y^ in-vector (generate (curry model beta y) y-sample-count))
                                 (summing (square (- (aref y i) y^)))))
                (vector-mean beta weight (/ 1 dev))))))

(defun likelihood (model y &key (offset 1))
  "Estimate P(y|model), assuming a normal distribution."
  (iter (for i from offset below (length y))
        (setf (fill-pointer y) i)
        (for samples = (generate (curry model y) 100))
        (multiply (normal-pdf (aref y i) (mean samples) (sd samples)))))

(defun posteriors (models priors y &key (offset 1))
  "Assign relative probabilities to each of models"
  (let* ((joint* (each (lambda (m p) (* p (likelihood m y :offset offset))) models priors))
         (marginal (sum joint*)))
    (e/ joint* marginal)))
