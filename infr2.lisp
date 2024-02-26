(require 'lisp-stat)
(require 'utils/misc)
(require 'iterate)
(require 'serapeum)

(defpackage :infr2
  (:shadowing-import-from :num-utils :product :left :right)
  (:shadowing-import-from :iter :next :generate :sum)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :lisp-stat :iter :plot :num-utils))

(in-package :infr2)

(defun wiener (n &key (var 1d0))
    (on (make-array n :element-type 'double-float)
      (iter (for x first 0d0 then (+ x (draw (r-normal 0d0 var))))
            (for i below n)
            (setf (aref this i) x))))

(defun simulate-wiener (dx/dt n &key (xi 0d0) (var 1d0))
  (on (make-array n :element-type 'double-float)
      (iter (for x first xi then (+ x (funcall dx/dt x) (draw (r-normal 0d0 var))))
            (for i below n)
            (setf (aref this i) x))))

(defun wiener-pdf (dx/dt y &key (var 1d0))
  (normal-pdf (elmt y -1) (+ (elmt y -2) (funcall dx/dt (elmt y -2))) var))

(defun posterior (likelihood prior marginal beta y &key (integrate #'romberg-quadrature))
  (/ (* (funcall likelihood beta y) (funcall prior beta))
     (funcall marginal likelihood prior y)))

(defun sample-marginal (n)
  (lambda (likelihood prior y)
    (mapn (draw prior) 100))
  )

(defun chain (model prior beta y* &optional (tail 1))
  (nlet chain ((prior prior) (i (1+ tail)))
    (print (funcall prior beta))
    (print (displace y* i))
    (if (> i (length y*))
        (funcall prior beta)
        (chain (lambda (beta) (posterior model prior beta (displace y* i))) (1+ i)))))

(defun my-chain (beta)
  (let* ((dx/dt (lambda (beta x) (* beta x)))
        (y* (simulate-wiener (curry dx/dt beta) 100))
        (model (lambda (beta y) (wiener-pdf (curry dx/dt beta) y)))
        (prior (lambda (beta) (normal-pdf beta 0d0 1d0))))
    (chain model prior beta y*)))

(my-chain 0.1)
