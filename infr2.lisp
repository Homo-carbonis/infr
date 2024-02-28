(require 'lisp-stat)
(require 'utils/misc)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr2
  (:shadowing-import-from :lisp-stat :product :next :sum :generate)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :lisp-stat :iter))

(in-package :infr2)

(defun wiener (n &key (var 1d0))
    (on (make-array n :element-type 'double-float)
      (iter (for x first 0d0 then (+ x (draw (r-normal 0d0 var))))
            (for i below n)
            (setf (aref this i) x))))

(defun simulate (step-x n &key (xi 0d0))
  (on (make-array n :element-type 'double-float :fill-pointer t)
      (iter (for x first xi then (funcall step-x x))
            (for i below n)
            (setf (aref this i) x))))

(defun step-wiener (dx/dt x &key (var 1d0))
  (+ x (funcall dx/dt x) (draw (r-normal 0d0 var))))

(defun wiener-pdf (dx/dt y &key (var 1d0))
  (normal-pdf (elmt y -1) (+ (elmt y -2) (funcall dx/dt (elmt y -2))) var))

(defun posterior (likelihood prior marginal beta y &key (integrate #'romberg-quadrature))
  (/ (* (pdf y (funcall likelihood beta)) (pdf beta prior))
     (funcall marginal likelihood prior y)))

; TODO: write mean reduction for iterate.

; Use nu.statistics:mean because lisp-stat:mean ignores weights if vector isn't just the right type.
(defun expected (model prior y &key (offset 0) (beta-samples 100) (y-samples 100) (mean #'nu.statistics:mean))
  "Find a Monte Carlo estimate of the parameter of a model.

   We draw a large number of samples from 'prior. Then for each sampled value
   of beta and each value of y generate samples from the model and compute the
   average squared difference from y. These are then used as weights to
   calculate an expected value for the parameter beta."
  (mean (iter (with beta* = (generate prior beta-samples))
              (for i from offset below (length y))
              (setf (fill-pointer y) i)
              (for weights = (iter (for beta in-vector beta*)
                             (for y^ = (mean (generate (curry model beta y) y-samples)))
                             (collect (/ 1 (square (- (aref y i) y^))))))
              (collect (funcall mean beta* :weights weights)))))

(defun model-likelihood (model y &key (offset 0))
  (mean (iter (for i from offset below (length y))
              (setf (fill-pointer y) i)
              (for samples = (generate (curry model y) 100))
              (collect (normal-pdf (aref y i) (mean samples) (sd samples))))))

(defun posteriors (models priors y &key (offset 0))
  (let* ((combined (map 'vector (lambda (m p) (* p (model-likelihood m y :offset offset))) models priors))
        (marginal (sum combined)))
    (map 'vector (lambda (c) (/ c marginal)) combined)))

(defun my-expected (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
        (y* (simulate (curry #'step-wiener (curry dx/dt beta)) 100))
        (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y))))
        (prior (curry #'draw (r-normal))))
    #++(vgplot:plot y*)
    (expected model prior y* :offset 1)))

(my-expected 0.1d0 )


(defun my-post (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
         (y* (simulate (curry #'step-wiener (curry dx/dt beta)) 100))
         (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y)))))
    #++(vgplot:plot y*)
    (posteriors (vector (curry model 0.0) (curry model 0.1) (curry model 0.2)) #(1/3 1/3 1/3) y* :offset 1)))

(my-post 0.1)
