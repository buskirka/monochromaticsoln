# monochromaticsoln
Code for resolving questions in Ramsey theory of the form

>Given k, f:N^j->R a linear function and a relation R (either =, <, or <=), 
>what is the greatest possible natural number m such that there exist colorings
>C: {1,...,m} -> {1,...,k} where if R(f(x_1,x_2,...,x_j),0) then it is not the 
>case that C(x_1)=C(x_2)=...=C(x_j) (a monochromatic solution)?

The `monochromaticsoln` function defines the core of the system; in its simplest usage
it takes the number of colors k and a vector specifying the function f, and returns a
matrix whose rows are the valid maximum-length k-colorings avoiding monochromatic
solutions to f(x_1,...,x_j)=0.
