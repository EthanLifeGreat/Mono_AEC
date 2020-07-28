# Mono_AEC

#### Description
Monophonic Acoustic Echo Cancellation (AEC) demo.

Using adaptive filter theory, find a filter w, according to input signals(reference signal x and desired signal d), that minimises (d - w * x), 

where d = x * w0 + r, r is Gaussian white noise and "*" denotes convolution operation


These AEC functions are copied and modified from Behrouz Farhang-Boroujeny (2013) Adaptive Filters:Theory and Applications 2nd Edition

#### Software Architecture
AECtest.m calls 5 AEC functions, namely:
1.    VSNLMS:
         Variable Step-size Nomalized LMS algorithm
2.    VSNLMSNt: 
         Variable Step-size Nomalized LMS-Newton algorithm
3.    VSAPLMS:
         Variable Step-size Affine Projection LMS algorithm
4.    VSNPFBLMS: 
         Variable Step-size Nomalized Partitioned Frequency-domain/Fast
         Block LMS algorithm
5.    SbLMS:
         Subband LMS algorithm

And, PFBfilter and SbFilter is the operation on x to create best-estimated d_hat using learned filter w.

#### Installation

1.  Clone the entire folder to Matlab's running path
2.  Run AECtest.m

#### Instructions

1.  Run AECtest.m and check the results
2.  Try changing algorithm options and check different algorithm's effect
3.  Try changing inputs (reference signal x and echo filter w0 ) and check the results 

