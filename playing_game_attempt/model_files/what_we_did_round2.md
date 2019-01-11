### Parameter Estimation

First run of differential evolution with ALL mRNA species included gave:
array([4.18050798e+00, 5.57396906e+00, 1.92921288e+00, 1.31558951e+00,
       6.48315203e+00, 4.58954214e+00, 4.24351002e-03])

Selecting for reasonable species (mRNA 1,3,4,6) gave: 
x: array([0.01556653, 9.959682  , 0.1056418 , 6.66957033, 0.08160472,
       4.25284957, 0.06687737])

Preliminary testing of both parameter sets seems to suggest the latter
is the better estimate.

Using mass spec data:
d_p = .037
a_p = 8.85

### Connection Predicting

First run of interaction estimation trying to predict connections to the orphans gene 2 and gene 7 gave:
Best connection for 2[(7, 2, -1), (8, 2, -1)]
Best connection for 7[(6, 7, -1), (7, 7, -1)]

### Parameter Picking

Vm1= 8.0
Vm2= 12.5
Vm3= 6.0
Vm4= 6
Vm5= 6
Vm6= 6
Vm7= 10
Vm8= 6

# Educated Guesses

Gene 2:
    SR: 7
    SR: 5
    SA: 8
Gene 8:
    DA: 6 + 8 
    DR 2 + 3
Gene 7:
    SR: 4
    DR: 2 + 4 
    DR: 5 + 7
Gene 6: 
(7, 6, -1)

[[(5, 7, -1), (7, 7, -1), (8, 2, -1), (3, 8, -1), (2, 8, 1)], 
[(7, 2, -1), (5, 7, -1), (3, 7, 1), (3, 8, -1), (2, 8, 1)], 
[(6, 8, 1), (8, 8, 1), (7, 2, -1), (5, 7, -1), (2, 7, 1)], 
[(5, 7, -1), (7, 7, -1), (6, 8, 1), (8, 8, 1), (5, 2, -1)], 
[(6, 8, 1), (8, 8, 1), (7, 7, -1), (5, 2, -1)], [(-1, -1, -1)]]
