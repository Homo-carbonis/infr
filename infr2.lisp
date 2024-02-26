(require 'lisp-stat)
(require 'utils/misc)
(require 'iterate)
(require 'serapeum)

(defpackage :infr2
  (:shadowing-import-from :num-utils :product :left :right)
  (:shadowing-import-from :iter :next :generate :sum)
  (:import-from :serapeum :nlet)
  (:use :cl :utils/misc :lisp-stat :iter :plot :num-utils)
  (:export :generate-noisy :posterior :posterior-2 :compose-normal-pdf :delta-normal-pdf))

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
  (normal-pdf (elmt y -1) (+ (elmt y -2) (funcall dx/dt (elt y -2))) var))

(defun posterior (likelihood prior beta y &key (integrate #'romberg-quadrature))
  (/ (* (funcall likelihood beta y) (funcall prior beta))
     (integrate (lambda (beta) (* (funcall likelihood beta y) (funcall prior beta))))))

(defun chain (model prior beta y* &optional (tail 1))
  (nlet ((chain ((prior prior) (i (+1 tail)))
                (if (> i (length y*))
                    (prior beta)
                    (chain (lambda (beta) (posterior model prior beta (displace y* i))) (1+ i)))))))




(chain (lambda (beta y) (wiener-pdf (curry dy/dx beta))))
