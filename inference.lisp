(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(defpackage :infr
  (:use :cl :utils/misc :lisp-stat :plot))

(in-package :infr)
(defparameter plot:*default-browser-command* :chrome-app-mode)

(defun generate-noisy (f variance n)
  (let ((y (make-array n :element-type 'double-float))
        (ti (make-array n :element-type 'double-float)))
    (loop for i from 0 below n
          do (setf (aref y i) (draw (r-normal (funcall f i) variance)))
          do (setf (aref ti i) (float  i 0.0d0)))
    (make-df '(:time :value) (list ti y))))

(defun f (beta ti) (exp (+ (* beta ti) (log 0.5))))

(defdf data (generate-noisy (curry #'f 0.1) 10 100))

(plot:plot
  (vega:defplot simple-line-plot
    `(:title "Noisy y data"
      :data (:values ,data)
      :mark :line
      :encoding (:x (:field :time
                     :type  :quantitative)
                 :y (:field :value
                     :type  :quantitative)))))

(defun likelihood (y ti beta sigma^2)
  (normal-pdf y (f beta ti) sigma^2))

(defun prior (beta &optional (mu 0) (gamma 1))
  (normal-pdf beta mu gamma))
