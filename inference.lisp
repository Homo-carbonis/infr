(require 'utils/misc)
(require 'alexandria)
(require 'plot/vega)
(require 'trivia)

(defpackage :infr
  (:shadowing-import-from :num-utils :product :left :right)
  (:use :cl :utils/misc :lisp-stat :plot :num-utils :trivia)
  (:export :generate-noisy :posterior :posterior-2 :compose-normal-pdf :delta-normal-pdf))

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

;TODO Redo with arrays. Write posterior-n with tail of y values.
(defun posterior-2 (pdf beta sigma prior y*)
  (let* ((deltas (deltas y*))
         (var (variance deltas))
         (mu (mean deltas)))
   (reduce-2 (lambda (prior y1 y2)
              (let ((marginal (normal-pdf (- y2 y1) mu var))
                    (likelihood (funcall pdf beta sigma y1 y2)))
                (format t "y1: ~a y2: ~a, prior: ~a~%" y1 y2 prior)
                (/ (* likelihood prior) marginal)))
            y*
            :initial-value prior)))

(defun compose-normal-pdf (f)
  (lambda (beta sigma x y)
    (normal-pdf y (funcall f beta x) sigma)))

(defun delta-normal-pdf (dx/dt)
  (lambda (beta sigma x1 x2)
    (normal-pdf x2 (+ x1 (funcall dx/dt beta x1)) sigma)))






