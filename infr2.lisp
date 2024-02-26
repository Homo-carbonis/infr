(require 'lisp-stat)
(require 'utils/misc)
(require 'iterate)
(require 'serapeum)
(require 'vgplot)

(defpackage :infr2
  (:shadowing-import-from :num-utils :product :left :right)
  (:shadowing-import-from :array-operations :generate)
  (:shadowing-import-from :iter :next :sum)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :lisp-stat :iter :plot :num-utils))

(in-package :infr2)

(defun wiener (n &key (var 1d0))
    (on (make-array n :element-type 'double-float)
      (iter (for x first 0d0 then (+ x (draw (r-normal 0d0 var))))
            (for i below n)
            (setf (aref this i) x))))

(defun simulate-wiener (dx/dt n &key (xi 0d0) (var 1d0))
  (on (make-array n :element-type 'double-float :fill-pointer t)
      (iter (for x first xi then (+ x (funcall dx/dt x) (draw (r-normal 0d0 var))))
            (for i below n)
            (setf (aref this i) x))))

(defun wiener-pdf (dx/dt y &key (var 1d0))
  (normal-pdf (elmt y -1) (+ (elmt y -2) (funcall dx/dt (elmt y -2))) var))

(defun posterior (likelihood prior marginal beta y &key (integrate #'romberg-quadrature))
  (/ (* (pdf y (funcall likelihood beta)) (pdf beta prior))
     (funcall marginal likelihood prior y)))

(defun posterior* (model prior beta y &optional (offset 0))
  (iter (with beta* = (generate (curry #'draw prior) 1000))
        (for i from offset below (length y))
        (setf (fill-pointer y) (1+ i))
        (for likelihood = (funcall model beta y))
        (for marginal = (mean (map-array beta* (lambda (beta) (funcall model beta y)))))
        (reducing (/ likelihood marginal) by #'* initial-value (pdf prior beta) into posterior)
        (format t "y: ~a, marginal: ~a, post: ~a~%" (elt y i) marginal posterior)
        (finally (return posterior))))




(defun my-posterior (beta beta^)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
        (y* (simulate-wiener (curry dx/dt beta) 100))
        (model (lambda (beta y) (wiener-pdf (curry dx/dt beta) y)))
        (prior (r-normal)))
    (vgplot:plot y*)
    (posterior* model prior beta^ y* 1)))

(my-posterior -0.1d0 -0.1d0)


