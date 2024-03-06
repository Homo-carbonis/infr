(require 'polisher)
(defpackage :infr-example
  (:import-from :misc-utils :elmt)
  (:use :cl :infr :lisp-stat :polisher :iterate))
(in-package :infr-example)
(polisher:activate-infix-syntax) 
(defun van-der-pol (mu sigma x)
  (let ( (h 1d-3)
        (x1 (elmt x -2))
        (x2 (elmt x -1)))
    (r-normal #i{(h* mu * x1 * x2**2 + (4 -2*h**2)*x2 - (h*mu +2) * x1) / (h*mu*x2**2 -h*mu+2)} sigma )))


(defparameter y (generate-markov-chain (curry #'van-der-pol 2d0 1d-9) 100000 :xi #(0.1d0 0.1d0)))

(vgplot:plot y)

(setf (fill-pointer y) 20000)

(exp (infr::log-posterior (lambda (beta y) (van-der-pol (first-elt beta) 1d-9 y))
                          (vector (r-normal 2d0 1d0))
                          #(3d0)
                          y
                          :sample-count 20))


(defparameter p-beta (iter (for mu from 0d0 to 5d0 by 0.5d0)
                          (collect (exp (infr::log-posterior
                                          #'van-der-pol
                                          (vector (r-normal 2d0 1d0) (r-uniform 0.5d-9 1.5d-9))
                                          #(1d0 1d-9)
                                          y
                                          :sample-count 50)))))
(vgplot:plot p-beta)
