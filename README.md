# monochromaticsoln
Code for resolving questions in Ramsey theory of the form

>Given a number of colors k, f:N^j->R a linear or affine function, 
>and a relation R (either =, <, or <=), 
>what is the greatest possible natural number m such that there exist colorings
>C: {1,...,m} -> {1,...,k} where if R(f(x_1,x_2,...,x_j),0) then it is not the 
>case that C(x_1)=C(x_2)=...=C(x_j) (a monochromatic solution)?

## Usage

The `monochromaticsoln` function defines the core of the system; in its simplest usage
it takes the number of colors k and a vector specifying the function f, and returns a
matrix whose rows are the valid maximum-length k-colorings avoiding monochromatic
solutions to f(x_1,...,x_j)=0. For example, to find all of the three-colorings to 
x_1+x_2=x_3, we may use
```
monochromaticsoln(3,[1,1,-1])
```
which specifies the number of colors to use (the 3) and the function ([1,1,-1], which 
is treated as x_1+x_2-x_3, which in the equation x_1+x_2-x_3=0 is equivalent to the 
original problem). The first two inputs to `monochromaticsoln` are *always* the number
of colors `k` (or a similar informative object) and the function specifier `EQN` 
as a row vector.

If we're working with a long-running problem, like finding 4-colorings of the equation
above, then we want more verbose output:
```
monochromaticsoln(4,[1,1,-1],'verbose')
```

To specify an affine function, as in finding MCS-avoiding colorings for x_1+x_2+1=x_3, 
we can use the `'constant'` flag:
```
monochromaticsoln(3,[1,1,-1],'Constant',1)
```
When Constant is defined, then the function f will be 
`EQN(1)*x_1+...+EQN(j)*x_j+L`.

To specify which relation we wish to use, for example to examine the inequality
x_1+x_2 < x_3 in 4 colors,
```
monochromaticsoln(4,[1,1,-1],'Relation','<')
```
Similarly, for x_1+x_2 <= x_3,
```
monochromaticsoln(4,[1,1,-1],'Relation','<=')
```
The opposite relations, `'>'` and `'>='` are not expected to be needed and thus are 
omitted from code, but may be obtained by negating `EQN` and `L` if needed for whatever 
reason.

Note that `monochromaticsoln` only finds solutions which introduce the colors in sequence;
that is, there is no 2 found before every 1, no 3 found before every 2, etc. To get every 
coloring in the broadest sense of the term, each possible permutation of the color set 
would need to be applied.
