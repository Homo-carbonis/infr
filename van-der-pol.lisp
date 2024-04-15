(defpackage :infr/van-der-pol
  (:use :cl :infr/utils :infr :lisp-stat :polisher :plot))

(in-package :infr/van-der-pol)
(polisher:activate-infix-syntax) 

; Finite difference approximatian of the Van der Pol oscillator with added Gaussian noise.
(defun van-der-pol (mu sigma x)
  (let ( (h 1d-3)
        (x1 (elmt x -2))
        (x2 (elmt x -1)))
    (r-normal #i{(h* mu * x1 * x2**2 + (4 -2*h**2)*x2 - (h*mu +2) * x1) / (h*mu*x2**2 -h*mu+2)} sigma )))

(defparameter f (lambda (beta y) (van-der-pol (elt beta 0) (elt beta 1) y)))

(defparameter y (generate-markov-chain (curry #'van-der-pol 2d0 1d-9) 100000 :xi #(0.1d0 0.1d0)))

(setf (fill-pointer y) (array-total-size y))

(defparameter prior (vector (r-normal 2d0 1d0) (r-uniform  0.9d-9 1.1d-9)))

(format t "log(P(y, beta =2)) ∝ ~a"
        (infr::log-likelihood f #(2d0 1d-9) y))

(format t "log(P(y | model)) ∝ ~a"
        (infr::log-marginal f
                            prior
                            y
                            :sample-count 50))

(format t "P(beta=2 | y) = ~a"
        (exp (infr::log-posterior f prior #(2d0 1d-9) y :sample-count 100)))
