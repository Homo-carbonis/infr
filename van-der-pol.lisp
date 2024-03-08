(require 'polisher)
(defpackage :infr-example
  (:import-from :misc-utils :elmt)
  (:import-from :misc-utils :elmt)
  (:shadowing-import-from :iterate :next)
  (:shadowing-import-from :lisp-stat :generate :sum)
  (:use :cl :infr :lisp-stat :polisher :iterate))
(in-package :infr-example)
(polisher:activate-infix-syntax) 
(defun van-der-pol (mu sigma x)
  (let ( (h 1d-3)
        (x1 (elmt x -2))
        (x2 (elmt x -1)))
    (r-normal #i{(h* mu * x1 * x2**2 + (4 -2*h**2)*x2 - (h*mu +2) * x1) / (h*mu*x2**2 -h*mu+2)} sigma )))

(defun van-der-pol (mu sigma x)
  (let ( (h 1d-3)
        (x1 (elmt x -2))
        (x2 (elmt x -1)))
    (r-normal #i{(h* mu * x1 * x2**2 + (4 -2*h**2)*x2 - (h*mu +2) * x1) / (h*mu*x2**2 -h*mu+2)} sigma )))

(defparameter f (lambda (beta y) (van-der-pol (first-elt beta) 1d-9 y)))

(defparameter y (generate-markov-chain (curry #'van-der-pol 2d0 1d-9) 100000 :xi #(0.1d0 0.1d0)))
(setf (fill-pointer y) (array-total-size y))

(vgplot:close-all-plots)
(vgplot:plot y)

(exp (infr::log-posterior f 
                          (vector (r-normal 2d0 1d0))
                          #(2d0)
                          y
                          :sample-count 2000))

(vgplot:plot (each #'exp (infr::log-posterior f 
                          (vector (r-normal 2d0 1d0))
                          #(2d0)
                          y
                          :sample-count 20 :seq t)))


(defparameter mu (linspace 0d0 5d0 100))
(defparameter p[mu] (each (lambda (mu)
                            (exp (infr::log-posterior
                                   f
                                   (vector (r-normal 2d0 1d0))
                                   (vector mu)
                                   y
                                   :sample-count 50)))
                          mu))
(vgplot:plot mu p[mu])
