(defpackage :example
  (:use :cl :infr :lisp-stat))

(in-package :example)

(defun dx/dt (beta x) (* beta x))

(defun generate-data (xi beta sigma n)
  (defdf data (generate-noisy (curry #'dx/dt beta) xi sigma n))
  (plot:plot
    (vega:defplot simple-line-plot
      `(:title "Noisy y data"
        :data (:values ,data)
        :mark :line
        :encoding (:x (:field :time
                       :type  :quantitative)
                   :y (:field :value
                       :type  :quantitative))))) )
 
(defun post (beta sigma y*)
 (posterior-2 (delta-normal-pdf #'dx/dt) beta sigma
              (normal-pdf beta 0.0 1.0) y*))

(generate-data 1d0 0.04d0 1d0 100)
(post 0.04 1d0 data:value)

(defun test-1 (test xi beta sigma n)
  (let* ((data (generate-noisy (curry #'dx/dt beta) xi sigma n)))
         (post test sigma data:value)))

(defun test (beta sigma)
  (let* ((beta* (make-array 2000 :element-type 'double-float))
         (p* (make-array 2000 :element-type 'double-float))
         (_ (loop for beta from (* -2 beta) below (* 2 beta) by (/ (* 4 beta) 2000)
                  for i from 0 below 2000
                  do (setf (aref beta* i) beta) 
                  do (setf (aref p* i) (post beta sigma data:value))))
         (beta-dist (make-df '(:beta :p) (list beta* p*))))
    #++(plot:plot
      (vega:defplot simple-line-plot
      `(:title "Data"
        :data (:values ,data)
        :mark :line
        :encoding (:x (:field :time
                       :type  :quantitative)
                   :y (:field :value
                       :type  :quantitative)))))
    (plot:plot
      (vega:defplot simple-line-plot
        `(:title "Beta distribution"
          :data (:values ,beta-dist)
          :mark :line
          :encoding (:x (:field :beta
                         :type  :quantitative)
                     :y (:field :p
                         :type  :quantitative)))))))
(test 0.04d0 1d0)
