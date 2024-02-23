(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(defpackage :infr
  (:use :cl :utils/misc :lisp-stat :plot :num-utils))

(in-package :infr)
(defparameter plot:*default-browser-command* :chrome-app-mode)


(defun generate-noisy (f variance n)
  (let ((y (make-array n :element-type 'double-float))
        (ti (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          do (setf (aref y i) (draw (r-normal (funcall f i) variance)))
          do (setf (aref ti i) (float  i 0.0d0)))
    (make-df '(:time :value) (list ti y))))

(defun example (xi beta ti) (exp (+ (* beta ti) (log xi))))

(defdf data (generate-noisy (curry #'example 0.5 0.1) 10 100))

(defun plot-data ()
 (plot:plot
  (vega:defplot simple-line-plot
    `(:title "Noisy y data"
      :data (:values ,data)
      :mark :line
      :encoding (:x (:field :time
                     :type  :quantitative)
                 :y (:field :value
                     :type  :quantitative))))))

(defun quadrature-2 (f interval-1 interval-2)
  (romberg-quadrature
    (lambda (x) (romberg-quadrature (curry f x) interval-2))
    interval-1))

(defun posterior (f y* ti* beta sigma^2)
  (let* ((pdf (lambda (y ti beta sigma^2) (normal-pdf y (funcall f beta ti) sigma^2)))
        (marginal (lambda (y ti)
                    (quadrature-2 (curry pdf y ti) (interval -1000 1000) (interval -1000 1000))))
        (likelihood (rcurry pdf sigma^2 beta))
        (prior (normal-pdf beta 0 1)))
    (/ (* (product (map 'vector likelihood y* ti*)) prior) (product (map 'vector marginal y* ti*)))))



(defun likelihood (y ti beta sigma^2)
  (normal-pdf y (f beta ti) sigma^2))

(defun prior (beta &optional (mu 0) (gamma 1))
  (normal-pdf beta mu gamma))

(defun marginal )
