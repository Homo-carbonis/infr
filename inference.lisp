(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(defpackage :infr
  (:use :cl :utils/misc :lisp-stat :plot :num-utils))

(in-package :infr)
(defparameter plot:*default-browser-command* :chrome-app-mode)


(defun generate-noisy (f variance n)
  (let ((ti (make-array n :element-type 'double-float))
        (y (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          do (setf (aref ti i) (float  i 0.0d0)) 
          do (setf (aref y i) (draw (r-normal (funcall f i) variance))))
    (make-df '(:time :value) (list ti y))))

(defun example (ti beta) (exp (+ (* beta ti) (log 0.5))))

(defdf data (generate-noisy (rcurry #'example 0.001) 1 100))

(defun plotdata ()
 (plot:plot
  (vega:defplot simple-line-plot
    `(:title "Noisy y data"
      :data (:values ,data)
      :mark :line
      :encoding (:x (:field :time
                     :type  :quantitative)
                 :y (:field :value
                     :type  :quantitative))))))

(defun quadrature-2 f interval-1 interval-2)

(defun quadrature-2 (f interval-1 interval-2)
  (romberg-quadrature
    (lambda (x) (romberg-quadrature (curry f x) interval-2 :epsilon 0.0001))
    interval-1 :epsilon 0.0001))

(defun posterior (pdf prior data* parameters)
  (reduce (lambda (prior data)
            (let* ((pdf (curry pdf data))
                   (marginal (quadrature-2 pdf (interval 0.11 0.2) (interval 0.1 0.2)))
                   (likelihood (apply pdf parameters)))
              (print marginal)
              (/ (* likelihood prior) marginal)))
          data*
          :initial-value prior))

(defun make-normal-over-f (f)
  (lambda (data sigma^2 &rest parameters)
    (let ((x (elt data 0))
          (y (elt data 1)))
    (normal-pdf y (apply f x parameters) sigma^2))))

(posterior (make-normal-over-f #'example)
           (normal-pdf 0.1 0 1) (rows data) '(0.001 10))

(defun i)
