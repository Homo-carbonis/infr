(require 'lisp-stat)
(require 'utils)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr2
  (:shadowing-import-from :lisp-stat :product :next :sum :generate :mean)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :iter-utils :lisp-stat :iter))

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

; Possibly we should weight the mean by magnitude of y?
(defun expected (model prior y &key (offset 0) (beta-sample-count 100) (y-sample-count 100) )
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

(defun model-likelihood (model y &key (offset 0) (chain-count 100))
  (generate (curry #'simulate model (length y) :xi (elt y 0) chain-count)
  (eexpt (e- y (simulate model (length y) :xi (elt y 0))) 2)
  

(defun posteriors (models priors y &key (offset 0))
  (let* ((joint* (each (lambda (m p) (* p (model-likelihood m y :offset offset))) models priors))
         (marginal (sum joint*)))
    (e/ joint* marginal)))

(defun my-expected (beta)
  (let* ((dx/dt (lambda (beta x) (+ (* (elt beta 0) x) (* (elt beta 1) (sqrt (abs x))))))
         (y* (simulate (curry #'step-wiener (curry dx/dt beta)) 100))
         (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y))))
         (prior (curry #'each #'draw (vector (r-normal) (r-normal 2d0)))))
    #++(vgplot:plot y*)
    (expected model prior y* :offset 1 :beta-sample-count 1000 :y-sample-count 100)))

(my-expected #(1d0 2d0))


(defun my-post (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
         (y* (simulate (curry #'step-wiener (curry dx/dt beta)) 100))
         (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y)))))
    #++(vgplot:plot y*)
    (posteriors (vector (curry model 0.09) (curry model 0.1) (curry model 0.11)) #(1/3 1/3 1/3) y* :offset 1)))

(my-post 0.101d0)
