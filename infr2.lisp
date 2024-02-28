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
  (:use :cl :utils/misc :lisp-stat :iter :num-utils))

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

(defun posterior* (model prior beta y &optional (offset 0))
  (iter (with beta* = (generate (curry #'draw prior) 1000))
        (for i from offset below (length y))
        (setf (fill-pointer y) (1+ i))
        (for likelihood = (funcall model beta y))
        (for marginal = (mean (map-array beta* (lambda (beta) (funcall model beta y)))))
        (reducing (/ likelihood marginal) by #'* initial-value (pdf prior beta) into posterior)
        (format t "y: ~a, marginal: ~a, post: ~a~%" (elt y i) marginal posterior)
        (finally (return posterior))))

(defun expected (model prior y &optional (offset 0))
  (mean (iter (with beta* = (generate (curry #'draw prior) 100))
              (for i from offset below (length y))
              (setf (fill-pointer y) i)
              (collect (iter (for beta in-vector beta*)
                             (for y^ = (mean (generate (curry model beta y) 100)))
                             (collect (/ 1 (square (- (aref y i) y^))) into weights)
                             (finally (return (nu.statistics:mean beta* :weights weights))))))))
;; lisp-stat:mean ignores weights if vector isn't just the right type

(defun posterior% (model prior y &optional (offset 0))
  (mean (iter (with beta* = (generate (curry #'draw prior) 100))
              (for i from offset below (length y))
              (setf (fill-pointer y) i)
              (for weights = (iter (for beta in-vector beta*)
                                   (for y^ = (mean (generate (curry model beta y) 100)))
                                   (collect (/ 1 (square (- (aref y i) y^))))))
              for beta^ = (mean beta* :weights weights)
              (map-array beta* (/ 1 (square (- (aref beta* i) beta^))))
              (collect beta^)
              )))

(defun my-expected (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
        (y* (simulate (curry #'step-wiener (curry dx/dt beta)) 100))
        (model (lambda (beta y) (step-wiener (curry dx/dt beta) (last-elt y))))
        (prior (r-normal)))
    #++(vgplot:plot y*)
    (expected model prior y* 1)))
(my-expected 0.1d0 )

(defun my-posterior (beta beta^)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
        (y* (simulate-wiener (curry dx/dt beta) 100))
        (model (lambda (beta y) (wiener-pdf (curry dx/dt beta) y)))
        (prior (r-normal)))
    (vgplot:plot y*)
    (posterior* model prior beta^ y* 1)))

(my-posterior -0.1d0 -0.1d0)


