(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(require 'drom/raster)

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

(defun quadrature-2 (f interval-1 interval-2 &optional (epsilon 0.001))
  (loop for i from (num-utils:left interval-1) below (num-utils:right interval-1) by epsilon
        sum (loop for j from (num-utils:left interval-2) below (num-utils:right interval-2) by epsilon
              sum (* epsilon (funcall f i j)))))

(defun posterior (pdf prior parameters)
  (let ((marginal (quadrature-2 pdf (interval 0.5 1.5) (interval -1.0 1.0)))
        (likelihood (apply pdf parameters)))
              (format t "prior: ~a, likelihood: ~a, marginal: ~a~%" prior likelihood marginal)
              (/ (* likelihood prior) marginal))
  )

(defun posterior (pdf prior data* parameters)
  (product a )
  (reduce (lambda (prior data)
            (let* ((pdf (curry pdf data))
                   (marginal (quadrature-2 pdf (interval 0.5 1.5) (interval -1.0 1.0)))
                   (likelihood (apply pdf parameters)))
              (format t "prior: ~a, likelihood: ~a, marginal: ~a~%" prior likelihood marginal)
              (/ (* likelihood prior) marginal)))
          data*
          :initial-value prior))

(defun posterior (pdf prior data* parameters)
  (reduce (lambda (prior data)
            (let* ((pdf (curry pdf data))
                   (marginal (quadrature-2 pdf (interval 0.5 1.5) (interval -1.0 1.0)))
                   (likelihood (apply pdf parameters)))
              (format t "prior: ~a, likelihood: ~a, marginal: ~a~%" prior likelihood marginal)
              (/ (* likelihood prior) marginal)))
          data*
          :initial-value prior))

(defun make-normal-over-f (f)
  (lambda (data sigma^2 &rest parameters)
    (let ((x (elt data 0))
          (y (elt data 1)))
    (normal-pdf y (apply f x parameters) sigma^2))))

(posterior (make-normal-over-f #'example)
           (normal-pdf 0.1 0.0 1.0) (rows data) '(1 0.001))
