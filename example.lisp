(defpackage :example
  (:use :cl :infr))
(defun example (beta ti) (exp (+ (* (first-elt beta) ti) (log 0.5))))

(defdf data (generate-noisy (curry #'example #(0.001)) 1 100))

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
 
(defun post (beta sigma^2)
 (posterior (make-normal-pdf-over-f #'example)
           (normal-pdf beta 0.0 1.0) (vector beta sigma^2) (rows data)))

(plot:plot
    (vega:defplot simple-line-plot
      `(:title "Noisy y data"
        :data (:values ,data)
        :mark :line
        :encoding (:x (:field :value
                       :type  :quantitative)
                   :y (:field :time
                       :type  :quantitative)))))

(let* ((beta* (make-array 2000 :element-type 'double-float))
      (p* (make-array 2000 :element-type 'double-float))
      (_ (loop for beta from -0.01d0 below 0.01d0 by 0.0001d0
              for i upfrom 0
              do (setf (aref beta* i) beta) 
              do (setf (aref p* i) (post beta 1)) ))
      (beta-dist (make-df '(:beta :p) (list beta* p*))))
  (plot:plot
    (vega:defplot simple-line-plot
      `(:title "Noisy y data"
        :data (:values ,beta-dist)
        :mark :line
        :encoding (:x (:field :beta
                       :type  :quantitative)
                   :y (:field :p
                       :type  :quantitative))))) )
