(defpackage :infr
  (:shadowing-import-from :lisp-stat :product :next :sum :generate :sd)
  (:shadowing-import-from :infr/utils :mean)
  (:use :cl :infr/utils :lisp-stat :iterate)
  (:export :generate-markov-chain :estimate-parameters :log-likelihood :log-marginal :log-posterior :posteriors))

(in-package :infr)

(defun generate-markov-chain (f n &key (xi #(0d0)))
  "Generate a realization of n elements of the Markov chain defined by f, beginning with xi."
  (iter (with offset = (length xi))
        (with x = (make-array offset
                              :initial-contents xi
                              :element-type 'double-float
                              :fill-pointer offset))
        (initially (adjust-array x n))
        (for i from offset below n)
        (vector-push (draw (funcall f x)) x)
        (finally (return x))))

; Possibly we should weight the mean by magnitude of y in some way?
(defun estimate-parameters (model prior y &key (offset 2) (sample-count 100))
  "Find a Monte Carlo estimate of the parameter of a model.  We draw a large
   number of samples from 'prior. Then for each sampled value of beta we find
   P(y|model,beta). Finally we take the mean of the samples using P(y|model,beta)
   as weights, to find an expected value for beta."
  (iter (with beta* = (generate (lambda () (draw* prior)) sample-count))
        (for i from offset below (length y))
        (setf (fill-pointer y) i) ; we pass y up to but not including element i to the model function.
        (mean
          (iter (for beta in-vector beta*)
                (vector-mean beta with-weight (pdf* (funcall model beta y) (aref y i)))))
        (finally (incf (fill-pointer y))))); Restore fill pointer before returning

(defun log-likelihood (model beta y &key (offset 2))
  "Find log P(y|model,beta)."
  (iter (for i from offset below (length y))
        (setf (fill-pointer y) i)
        (summing (log-pdf (funcall model beta y) (aref y i)))
        (finally (incf (fill-pointer y)))))

(defun log-marginal (model prior y &key (offset 2) (sample-count 100))
  "Estimate log P(y|model). We take a random sample of beta and use it to estimate the integral of P(y|beta)P(beta) over beta."
  ;; We factor out the first term, likelihood0, to compute the arithmetic mean
  ;; from logs, weighted by 1/prior.
  (iter (with beta* = (generate (lambda () (draw* prior)) sample-count))
          (with likelihood0 = (log-likelihood model (elt beta* 0) y :offset offset))
          (for beta in-vector beta*)
          (for likelihood = (log-likelihood model beta y :offset offset))
          (for p = (pdf* prior beta))
          (for term first 1d0 then (exp (- likelihood likelihood0)))
          (summing term into sum)
          (summing (/ 1 p) into weight)
          (finally (return (+ (- (log weight))
                              (sample-hypervolume beta*)
                              likelihood0
                              (log sum))))))

(defun log-posterior (model prior beta y &key (offset 2) (sample-count 100))
  "Calculate log p(beta| model, y)"
  (let ((likelihood (log-likelihood model beta y :offset offset))
        (prior (log-pdf* prior beta))
        (marginal (log-marginal model prior y :offset offset
                                :sample-count sample-count)))
    (+ likelihood prior 
       (- marginal))))

(defun pdf* (rv* x*)
  (product (each (lambda (rv x) (pdf rv x)) rv* x*)))

(defun log-pdf* (rv* x*)
  (product (each (lambda (rv x) (log-pdf rv x)) rv* x*)))

(defun draw* (distribution)
  (each #'draw distribution))

(defun min* (a b)
  (vectorize (a b) (min a b)))

(defun max* (a b)
  (vectorize (a b) (max a b)))

(defun sample-hypervolume (sequence)
  "Find the minimum volume containing every sample in a sequence of n-vectors"
  (let ((min (reduce #'min* sequence))
        (max (reduce #'max* sequence)))
    (product (vectorize (max min) (- max min)))))
