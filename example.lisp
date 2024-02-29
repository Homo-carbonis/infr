(defpackage :infr-example
  (:use :cl :infr))
(in-package :infr-example)
(defun my-expected (beta)
  (let* ((dx/dt (lambda (beta x) (+ (* (elt beta 0) x) (* (elt beta 1) (sqrt (abs x))))))
         (y* (generate-markov (curry #'step-wiener (curry dx/dt beta)) 100))
         (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y))))
         (prior (curry #'each #'draw (vector (r-normal) (r-normal 2d0)))))
    #++(vgplot:plot y*)
    (estimate-parameters model prior y* :offset 1 :beta-sample-count 1000 :y-sample-count 100)))

(my-expected #(1d0 2d0))


(defun my-post (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
         (y* (generate-markov (curry #'step-wiener (curry dx/dt beta)) 100))
         (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y)))))
    #++(vgplot:plot y*)
    (posteriors (vector (curry model 0.09) (curry model 0.1) (curry model 0.11)) #(1/3 1/3 1/3) y* :offset 1)))

(my-post 0.11d0)
