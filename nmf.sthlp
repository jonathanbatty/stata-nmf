{smcl}
{* *! version 1.0  16sep2022}{...}
{viewerjumpto "Syntax" "nmf##syntax"}{...}
{viewerjumpto "Description" "nmf##description"}{...}
{viewerjumpto "Options" "nmf##options"}{...}
{viewerjumpto "Remarks" "nmf##remarks"}{...}
{viewerjumpto "References" "nmf##references"}{...}
{viewerjumpto "Examples" "nmf##examples"}{...}
{title:Title}

{phang}
{bf:nmf} {hline 2} matrix decomposition using non-negative matrix factorization (NMF).


{marker syntax}{...}
{title:Syntax}

Perform matrix decomposition using NMF:

{p 8 17 2}
{cmdab:nmf}
[{varlist}]
{cmd:,} {bf: k(#)} {bf: epoch(#)} [{it:options}]

{synoptset 20 tabbed}{...}
{synopthdr}
{synoptline}
{syntab:Main}
{p2coldent:* {opt k(#)}}perform matrix decomposition to give non-negative matrices of rank {bf: k(#)}.{p_end}
{p2coldent:* {opt epoch(#)}}maximum number of iterations over which to optimize the decomposed matrices W and H to minimise the error function.{p_end}
{synopt:{opt initial(string)}}used to set initialisation method. Options include: randomu, randomn, nndsvd, nndsvda and nndsvdar.{p_end}
{synopt:{opt loss(string)}}defines the loss function used to calculate the error function for {c |}{c |}A - WH{c |}{c |}. Options include: eu, is and kl. Default is {cmd: loss(eu)}{p_end}
{synopt:{opt stop(#)}}sets the early stopping threshold for convergence; if {cmd: stop(0)} is set, optimization will continue for the set number of epochs. If ((previous error - current error) / error at initiation) < stop tolerance), convergence has occured and optimization terminates. Default is {cmd: stop(1.0e-4)}{p_end}
{synopt:{opt nograph}}plots graph of epoch vs. the chosen loss function at the end of the decomposition.{p_end}
{synoptline}
{p2colreset}{...}
{p 4 6 2}
* {opt k(#)} is required.


{p 4 6 2}
{cmd:by} is not allowed. {cmd:fweight}s are not allowed.



{marker description}{...}
{title:Description}

{pstd}
The aim of {cmd:nmf} is to decompose matrix A into a low-rank matrix approximation such that 
X ≈ WH, where A, W and H are made up of non-negative, real numbers. Missing values are 
permitted. If X is an n x m matrix, the dimensions of W and H will be n x k and k x m, 
respectively, whereby k represents the rank of the decomposition. In many cases, a good 
approximatation for A may be achieved with k << rank(A). 

{pstd}
NMF is a NP-hard problem. As such, multiplicative updates of matrices W and H are iteratively  
performed in order to minimise the generalized error function, {c |}{c |}A - WH{c |}{c |}.


{marker options}{...}
{title:Options}

{dlgtab:Main}

{phang}
{opt k(#)} is required and specifies the rank for the decomposition.

{phang}
{opt epoch(#)} is required and sets the maximum number of epochs (iterations) over which the decomposition is optimized.

{dlgtab:Options}

{marker initial()}{...}
{phang}
{opt initial(option)}
indicates how matrices {it: W} and {it: H} are initialized.  The
available {it:option}s are:

{phang2}
{opt randomu},
the default, specifies that .

{phang2}
{opt randomn},
specifies that .

{phang2}
{opt nndsvd},
specifies that .

{phang2}
{opt nndsvda},
specifies that 

{phang2}
{opt nndsvdar},
specifies that .

{marker loss()}{...}
{phang}
{opt loss(option)}
indicates which loss function will be minimized during the optimization of the decomposition.  The
available {it:option}s are:

{phang2}
{opt is},
the default, specifies that .

{phang2}
{opt eu},
specifies that .

{phang2}
{opt kl},
the default, specifies that .

{marker stop()}{...}
{phang}
{opt stop(#)}
indicates the tolerance for early stopping.  The
available {it:option}s are:

{phang2}
{opt 0},
specifies that .

{phang2}
{opt 1.0e-4},
the default, specifies that .

{marker nograph}{...}
{phang}
{opt nograph}
indicates the tolerance for early stopping.  The
available {it:option}s are:



{marker remarks}{...}
{title:Remarks}

{pstd}
For detailed information on the whatever statistic, see
{manlink R Intro}.

{marker references}{...}
{title:References}

{pstd}
[1] Paatero, P. and Tapperm U. Positive matrix factorization: A non-negative factor model with optimal utilization of error estimates of data values. Environmetrics, 5: pp. 111-126 (1994).
{p_end}

{pstd}
[2] Lee, D. and Seung, H. Learning the parts of objects by non-negative matrix factorization. Nature 401, pp. 788–791 (1999).
{p_end}

{pstd}
[3] Lee, D. and Seung, H. Algorithms for Non-negative Matrix Factorization. Advances in Neural Information Processing Systems 13: Proceedings of the 2000 Conference, MIT Press. pp. 556–562 (2000).
{p_end}

{marker examples}{...}
{title:Examples}

{phang}{cmd:. nmf p*, k(5) epoch(100)}{p_end}

{phang}{cmd:. nmf k*,	k(15) epoch(100) initial(randomu) stop(1.0e-4) loss(kl) nograph	}{p_end}

