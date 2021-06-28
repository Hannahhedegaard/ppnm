Examination project: Practical Programming and Numerical Methods 2021.
Hannah Hedegaard Nielsen, 201704744. Email: hannah@chem.au.dk

My examination project is: 44 modulo 22 = 0 -> Akima sub-spline implementation.

Note: all references to figures, tables, and equations are to the book: 
[1] "Yet Another Introduction to Numerical Methods" version May 12, 2021, D.V. Fedorov

--------------------------------------------------------------------------------------
Introduction

If a function only is represented by a discrete set of data-points, a smooth
function can be constructed, which passes exactly through the data-points and will
approximate the function between these points. This is called interpolation.
Interpolation is thus a great tool for estimating values between known data-points.
The interpolating function can vary - often is spline interpolation used, where a
piecewise polynomial, called a spline, is used as the interpolating function. The
most used spline is the cubic spline, which is made of third order polynomials.

In this examination project, the Akima sub-spline is implemented. Sub-splines also
uses piecewise polynomials as the interpolating function, but they do not have 
conditions for maximal differentiability of the spline.
When the data contains an outlier point, the splines often display unpleasant 
wiggles. Sub-splines try to minimize these wiggles.

The Akima sub-spline is a piecewise cubic polynomial, and is therefore similar to
the cubic spline. The equation is shown in (1.30). The difference between these
two is, that the Akima sub-spline does not have a requirement of maximal 
differentiability, which means the second derivative is discontinued. 
The advantage of this sub-spline is that it reduces wiggling seen for the cubic
spline. [1]

--------------------------------------------------------------------------------------
Implementation

The Akima sub-spline is implemented in the file functions.c. Three different
functions are implemented, one to allocate vectors with the sub-spline coefficients,
one to evaluate the value of the spline in a given point, and one to free memory.

The implementation is based on the equations given in chapter 1.3.4, together with
the implemented functions in Homework 1 and Table 1.5.

The different steps in the implementation is commented in the functions.c file.

--------------------------------------------------------------------------------------
Testing and further comments

The implemented Akima sub-spline is tested on two different data-sets.
One is to make a replication of Figure 1.2 [1], and the other is to show how the
sub-spline interpolates when there is one outlier in the data.
In both examples, the Akima sub-spline is compared to the cubic spline. This spline
was implemented in Homework 1, and the functions are thus taken therefrom.

The out.txt file contains comments on these two examples and shows the results from
the implementation - together with the two figures, Figure1 and Figure2.

--------------------------------------------------------------------------------------