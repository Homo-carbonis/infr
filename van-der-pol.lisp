(defpackage :infr/van-der-pol
  (:shadowing-import-from :iterate 
  (:shadowing-import-from :lisp-stat :generate :sum :mean)
  (:use :cl :infr/utils :infr :lisp-stat :polisher :iterate :plot)))
(in-package :infr/van-der-pol)
(polisher:activate-infix-syntax) 

; Finite difference approximatian of the Van der Pol oscillator with added Gaussian noise.
(defun van-der-pol (mu sigma x)
  (let ( (h 1d-3)
        (x1 (elmt x -2))
        (x2 (elmt x -1)))
    (r-normal #i{(h* mu * x1 * x2**2 + (4 -2*h**2)*x2 - (h*mu +2) * x1) / (h*mu*x2**2 -h*mu+2)} sigma )))

; We treat sigma as a constant for the sake of a simple 1d example.
(defparameter f (lambda (beta y) (van-der-pol beta 1d-9 y)))

(defparameter y (generate-markov-chain (curry #'van-der-pol 2d0 1d-9) 100000 :xi #(0.1d0 0.1d0)))
(setf (fill-pointer y) (array-total-size y))
(defdf y-data (make-df '(:y :time) (list y (linspace 0 (1- (length y)) (length y)))))

(plot:plot
  (vega:defplot van-der-pol
  `(:title "Van Der Pol Oscillator"
    :description ""
    :data (:values ,y-data)
    :mark :line
    :encoding (:y (:field :y :type :quantitative)
               :x (:field :time :type :quantitative)))))

(defdf marginal-data (infr::log-marginal f 
                    (r-normal 2d0 1d0)
                    y
                    :sample-count 50 :plot t))

(plot:plot
  (vega:defplot marginal
  `(:title "Particle estimation of marginal"
    :description "Particle estimation of the marginal probability of y. We sample values of beta from the prior distribution and take a weighted mean of the marginals to estimate the integral."
    :data (:values ,marginal-data)
    :mark :bar
    :encoding (:x (:field :beta :type :quantitative)
               :y (:field :posterior :type :quantitative :stack "none" :scale (:zero :false))))))

(format t "log(P(y, beta =2)) ∝ ~a"
        (infr::log-likelihood f 2d0 y))

(format t "log(P(y | model)) ∝ ~a" (last-elt marginal-data:estimate))

(format t "P(beta=2 | y) = ~a"
        (exp (infr::log-posterior f (r-normal 2d0 1d0) 2d0 y :sample-count 100)))
