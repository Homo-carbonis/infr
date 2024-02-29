(require 'lisp-stat)
(require 'utils)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr
  (:shadowing-import-from :lisp-stat :product :next :sum :generate :mean)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :iter-utils :lisp-stat :iter)
  (:export :generate-wiener :generate-markov :step-wiener :estimate-parameters :likelihood :posteriors))

(in-package :infr)

(defun generate-markov (f n &key (xi 0d0))
  (on (make-array n :element-type 'double-float :fill-pointer t)
      (iter (for x first xi then (funcall f x))
            (for i below n)
            (setf (aref this i) x))))

(defun step-wiener (dx/dt x &key (var 1d0))
  (+ x (funcall dx/dt x) (draw (r-normal 0d0 var))))

; Possibly we should weight the mean by magnitude of y in some way?
(defun estimate-parameters (model prior y &key (offset 0) (beta-sample-count 100) (y-sample-count 100) )
  "Find a Monte Carlo estimate of the parameter of a model.
   We draw a large number of samples from 'prior. Then for each sampled value
   of beta and each value of y generate samples from the model and compute the
   average squared difference from y. These are then used as weights to
   calculate an expected value for the parameter beta."
  (iter (with beta* = (generate prior beta-sample-count))
        (for i from offset below (length y)) 
        (setf (fill-pointer y) i) ; we pass y up to but not including element i to the model function.
        (vector-mean
          (iter (for beta in-vector beta*)
                (for dev = (iter (for y^ in-vector (generate (curry model beta y) y-sample-count))
                                 (summing (square (- (aref y i) y^)))))
                (vector-mean beta weight (/ 1 dev))))))

(defun likelihood (model y &key (offset 0))
  "Estimate P(y|model), assuming a normal distribution."
  (iter (for i from offset below (length y))
        (setf (fill-pointer y) i)
        (for samples = (generate (curry model y) 100))
        (multiply (normal-pdf (aref y i) (mean samples) (sd samples)))))

(defun posteriors (models priors y &key (offset 0))
  "Assign relative probabilities to each of models"
  (let* ((joint* (each (lambda (m p) (* p (likelihood m y :offset offset))) models priors))
         (marginal (sum joint*)))
    (e/ joint* marginal)))
