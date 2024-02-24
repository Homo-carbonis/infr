(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(require 'drom/raster)

(defpackage :infr
  (:shadowing-import-from :num-utils :product :left :right)
  (:use :cl :utils/misc :lisp-stat :plot :num-utils)
  (:export :generate-noisy :posterior :make-normal-pdf-over-f))

(in-package :infr)

(defun generate-noisy (dx/dt xi variance n)
  (let ((ti* (make-array n :element-type 'double-float))
        (y* (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          for y = xi then (+ y (draw (r-normal (funcall dx/dt y) variance)))
          do (setf (aref ti* i) (float  i 0.0d0)) 
          do (setf (aref y* i) y))
    (make-df '(:time :value) (list ti* y*))))

(defun posterior (pdf prior parameters data*)
  (reduce (lambda (prior data)
            (let* ((marginal 1)
                   (likelihood (funcall pdf parameters data)))
              (/ (* likelihood prior) marginal)))
          data*
          :initial-value prior))

(defun make-normal-pdf-over-f (f)
  (lambda (parameters data)
    (let ((x (elt data 0))
          (y (elt data 1))
          (sigma^2 (last-elt parameters))
          (parameters (subseq parameters 0 (1- (length parameters))))) 
    (normal-pdf y (funcall f parameters x) sigma^2))))







