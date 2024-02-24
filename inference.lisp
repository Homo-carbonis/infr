(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(require 'drom/raster)

(defpackage :infr
  (:use :cl :utils/misc :lisp-stat :plot :num-utils :drom/raster :drom/interval))

(in-package :infr)
(defparameter plot:*default-browser-command* :chrome-app-mode)


(defun generate-noisy (f variance n)
  (let ((ti (make-array n :element-type 'double-float))
        (y (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          do (setf (aref ti i) (float  i 0.0d0)) 
          do (setf (aref y i) (draw (r-normal (funcall f i) variance))))
    (make-df '(:time :value) (list ti y))))


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







