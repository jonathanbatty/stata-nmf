capture prog drop nmf
*! version 0.01 05 Sep 2022
program define nmf, rclass
    version 17
 
    syntax varlist(numeric), k(integer) iter(integer) [initial(string) beta(numlist integer max = 1) stop(numlist max = 1) nograph noframes]
    
    // Aim:
    // The aim of NMF is to find W (u x k) and H (k x v) such that A ~ WH, where all 3 matrices contain only non-negative values
    // A good appoximation may be achieved with k << rank(A)

    // Parameters:
    // Initialisation options, specified using the init() parameter, include:
    //    random [*default] - random initialisation of matrix in range [0, 1] 
    //    nndsvd - Nonnegative Double Singular Value Decomposition [Boutsidis2007] - suited to sparse factors
    //    nndsvda - NNSVD with zero elements replaced with the input matrix average (not recommended for sparse data)
    //    nndsvdar - NNSVD with zero elements replaced with a random value (nor recommended for sparse data)
    // Using the random option, multiple runs will be required to identify the factorization that achieves the lowest approximation error. Other options are deterministic and only need to be ran once.
    // Note: The multiplicative update ('mu') solver cannot update zeros present in the initialization, and so leads to poorer results when used jointly with nndsvd. Co-ordinate descent mitigates this

    // Beta divergence options, specified using the beta() parameter, include:
    //    0 - Itakura-Saito divergence 
    //    1 - Generalized Kullback-Leibler divergence
    //    2 [*default] - Frobenius (Euclidean) norm

    // Early stopping options, specified using the stop() parameter, include:
    //    0 - early stopping off; NMF continues for numer of iterations specified in iter
    //    float value, e.g. 1.0e-4 [*default] - if ((previous error - current error) / error at initiation) < stop tolerance) then further iterations are not performed

    // nograph suppresses production of the graph of beta divergence at each iteration
    // noframes does not save frames containing output matrices: W, H and norms 

    // Uses multiplicative update method of Lee and Seung in 1999: 
    // https://proceedings.neurips.cc/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf

    // To assess the ideal value of k:
    //      Plot k vs residual sum of squares (i.e. Frobenius norm) and look for the value of r that shows an inflection point
    //      Plot k vs cophenetic correlation coefficient
    // 

    // To do: 
    // implement solver() - Alternate Least Square (ALS) approach of co-ordinate descent,  Kim H, Park H: Sparse non-negative matrix factorizations via alternating non-negativity-constrained least squares for microarray data analysis. Bioinformatics (Oxford, England) 2007, 23:1495-502 [http://www.ncbi.nlm.nih.gov/pubmed/17483501]
    // Implement cophenetic correlation coefficient as a measure of stability of the clusters -  Brunet JP, Tamayo P, Golub TR, Mesirov JP: Metagenes and molecular pattern discovery using matrix factorization. Proceedings of the National Academy of Sciences of the United States of America 2004, 101:4164-9 [http://www.ncbi.nlm.nih.gov/pubmed/15016911].
    // Consider scaling random values -? by sqrt(average of input matrix / length of input matrix) as per skl?
    // Consensus clustering: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC384712/
    //    1. After each run of NMF, allocate each patient to a 'cluster' based on the highest value in the result matrix 
    //    2. Plot a connectivity matrix, C, of M x M (patients x patients) with 1 = in same cluster and 0 = not in same cluster. Save this matrix.
    //    3. Average the connnectivity matrices over all runs (number of runs is selected by continuing until C̄ appears to stabilize; typically 20–100 runs suffice). Entries 
    //       of C are in the range [0, 1]. This dispersion between 0 and 1 therefore measures the reproducibility of the class assignments with respect to random initial conditions
    //       This wll eb a symmerical matrix, diagonals will be 1.
    //    4. Use average linkage HC to reorder the patients and thus the rows and columns of C
    //    5. Evaluate the stability of clustering associated with a given rank k
    //    6. Calculate the cophenetic correlation coefficient ρk(C̄) which indicates the dispersion of the consensus matrix C̄
    //         (a) ρk is computed as the Pearson correlation of two distance matrices: the first, I-C̄, is the distance between samples induced by the consensus matrix, 
    //             and the second is the distance between samples induced by the linkage used in the reordering of C̄
    //       In a perfect consensus matrix (all entries = 0 or 1), the cophenetic correlation coefficient equals 1. When the entries are scattered between 0 and 1, the cophenetic correlation coefficient is <1.
    //       We observe how ρk changes as k increases. We select values of k where the magnitude of the cophenetic correlation coefficient begins to fall

    // Consider adding a regularisation function defined to enforce, e.g. smoothness, sparsity on W and H
    // Consider implementing a version for large matrices (which uses st_view instead of st_data and minimises copy operations, uses cross() for matmul, etc)
    // Test and benchmark
    // 

    // If a value for initialisation method is not passed, default option is is random initialisation 
    if "`initial'" == "" local initial "random"

    // If a value specifying the form of beta divergence normalisation is not passed, default option is is Frobenius normalisation (2) 
    if "`beta'" == "" local beta = 2

    // If a value specifying the stopping condition is not passed, default option is is 1.0e-4 (i.e. 0.01%)
    if "`stop'" == "" local stop = 1.0e-4

    // Run NMF
    mata: nmf("`varlist'", `k', `iter', "`initial'", `beta', `stop')


    // Store using stata matrices
    matrix W = r(W)
    matrix H = r(H)
    matrix norms = r(norms)

    // Creates frames containing output matrices: W, H and norms
    if "`frames'" != "noframes" {
        
        display "Creating frames W, F and norms to hold results."

        // Create empty new frames to hold matrices
        foreach outputFrame in W H norms {

            // If a frame with the given name exists - drops and recreates it
            capture confirm frame `outputFrame'
            if !_rc {
                frame drop `outputFrame'
                frame create `outputFrame'
            }
            else {
                frame create `outputFrame'
            }

            // Populate each frame with the returned matrix object
            frame `outputFrame' {
                quietly svmat `outputFrame'
            }
        }

        // Add iteration identifier to norms frame
        frame norms {
            gen iteration = _n
            order iteration
        }
    }

    // Plots the normalisation values calculated after each iteration
    if "`graph'" != "nograph" {

        display "Plotting graph of loss function..."

        if `beta' == 0 local betaString "Itakura-Saito Divergence" 
        if `beta' == 1 local betaString "Generalized Kullback-Leibler Divergence"
        if `beta' == 2 local betaString "Frobenius (Euclidean) Normalization"

        local plotIterations = colsof(norms)

        frame norms: graph twoway line norm iteration, ///
        title("Loss Function") ///
        xtitle("Iteration") xlabel(, labsize(*0.75) grid glcolor(gs15)) ///
        ytitle("`betaString'") yscale(range(0 .)) ylabel(#5, ang(h) labsize(*0.75) grid glcolor(gs15)) ///
        graphregion(color(white))
    }

    display "Returning final matrices in r(W), r(H) and r(norms)"
    // Return final results
    return matrix W W
    return matrix H H
    return matrix norms norms

end

version 17
mata:

void nmf(string scalar varlist, 
         real scalar k, 
         real scalar iter, 
         string scalar initial,
         real scalar beta,
         real scalar stop)
{
    // Construct mata matrix from input varlist
    A = st_data(., varlist)

    // Perform error checking:
    // 1. Ensure that input matrix is non-negative (i.e. no negative values)
    if (min(A) < 0) {
        _error("The input matrix must not contain negative values.")
    }

    // 2. Ensure that there are no missing values in the input matrix
    if (hasmissing(A) == 1) {
        _error("The input matrix must not contain missing values.")
    }

    // 3. Ensure that specified rank, k,  is valid (i.e. 2 < k < rank(A))
    if (k < 2 | k >= cols(A) | k >= rows(A)) {
        _error("The rank of the output matrix must be at least 2, but less than the number of rows and columns of the input matrix, A.")
    }

    // Declaration of output matrices
    real matrix W, H

    // Initialisation of matrix
    if (initial == "random") {
        randomInit(A, k, W, H)
    }
    else if (initial == "nndsvd" | initial == "nndsvda" | initial == "nndsvdar") {
        nnsvdInit(A, k, W, H, initial)
    }
    else {
        printf("Error - an invalid initialisation type has been set.")
    }

    // Declare value of e
    e = 1.0e-10
    
    // Create a column matrix (rows = iterations, columns = 1) to store normalisation results from each iteration
    real matrix norms 
    norms = J(1, 1, .)

    // Perform multiplicative updating given number of iterations
    printf("Factorizing matix using multiplicative update method.\n\n")
    multiplicativeUpdating(A, k, iter, W, H, e, norms, beta, stop)

    // Return results object to stata
    st_matrix("r(W)", W)
    st_matrix("r(H)", H)
    st_matrix("r(norms)", norms)
}


void randomInit(real matrix A, 
                real scalar k, 
                real matrix W, 
                real matrix H)
{
    // Generate random values for W in the range [1 - 2]
    W = runiform(rows(A), k, 1, 2)

    // Generate random values for W in the range [1 - 2]
    H = runiform(k, cols(A), 1, 2)
}

void nnsvdInit(real matrix A, 
               real scalar k, 
               real matrix W, 
               real matrix H, 
               string scalar initial)
{
    // Declare empty matrices to hold results of SVD
    real matrix U
    real matrix S
    real matrix V

    // Perform SVD and transpose resulting S matrix
    svd(A, U, S, V)
    S = S'

    // Set up empty matrices of the correct size for W and H
    W = J(rows(A), k, 0)
    H = J(k, cols(A), 0)
    
    // Get first column of W values based on SVD results
    W[., 1] = sqrt(S[1, 1]) :* abs(U[., 1])

    // Get first row of H values based on SVD results
    H[1, .] = sqrt(S[1, 1]) :* abs(V[., 1]')

    for (i = 2; i <= k; i++) {

        ui = U[., i]
        vi = V[., i]

        // Divide into positive-only and negative-only matrices
        ui_pos = (ui :>= 0) :* ui
        ui_neg = (ui :< 0) :* -ui
        vi_pos = (vi :>= 0) :* vi
        vi_neg = (vi :< 0) :* -vi

        // Calculate 2-norm of each of the positive and negative columns
        ui_pos_norm = norm(ui_pos, 2)
        ui_neg_norm = norm(ui_neg, 2)
        vi_pos_norm = norm(vi_pos, 2)
        vi_neg_norm = norm(vi_neg, 2)

        // Multiply the positive and negative norms to get overall values
        norm_pos = ui_pos_norm * vi_pos_norm
        norm_neg = ui_neg_norm * vi_neg_norm

        // Check which is larger, norm_pos or norm_neg
        if (norm_pos >= norm_neg) {
            _ui = ui_pos :/ ui_pos_norm
            _vi = vi_pos :/ vi_pos_norm
            sigma = norm_pos
        }
        else {
            _ui = ui_neg :/ ui_neg_norm; 
            _vi = vi_neg :/ vi_neg_norm;
            sigma = norm_neg
        }
        
        // Update rows of W and H
        W[., i] = sqrt(S[1, i] * sigma) * _ui
        H[i, .] = sqrt(S[1, i] * sigma) * _vi'
    }

    // Calculate properties of input matrix for alternative methods below
    averageOfInputMatrix = length(A) / sum(A)

    // Replaces zero values in initialised matrices with average value of input matrix, A
    if (initial == "nndsvda") {
        
        // Update zeros with average values of input matrix
        W = editvalue(W, 0, averageOfInputMatrix)
        H = editvalue(H, 0, averageOfInputMatrix)

    }
    // Replace zeros in initialised matrices with random value in the space [0 : average/100]
    else if (initial == "nndsvdar") {

        // Update zeroes in W with averge vales, scaled
        W = editvalue(W, 0, averageOfInputMatrix * runiform(1, 1) / 100)
        H = editvalue(H, 0, averageOfInputMatrix * runiform(1, 1) / 100)
    }
    // Replace zeros in initialised matrices with the smallest possible positive value
    else {
        W = editvalue(W, 0, epsilon(1))
        H = editvalue(H, 0, epsilon(1))
        
    }
}

scalar betaDivergence(real matrix A,
                      real matrix W,
                      real matrix H,
                      real scalar beta)
{
    if (beta == 0) {
        // 0 = Itakura-Saito divergence (only if no zero/missing values)
        div = A :/ (W*H)
        divergence = sum(div) - (rows(A) * cols(A)) - sum(log(div))
    }
    else if (beta == 1) {
        // 1 = Generalized Kullback-Leibler divergence
        WH = W*H
        
        // Unravel A and WH to single row vectors
        A_data = rowshape(A, 1)
        WH_data = rowshape(WH, 1)

        // Replace values of WH that are 0 to a very small value to prevent div by 0 errors
        A_data = (A_data :>= epsilon(1)) :* A_data
        WH_data = (A_data :>= epsilon(1)) :* WH_data

        WH_data = editvalue(WH_data, 0, epsilon(1))
        
        sum_WH = sum(WH)
        div = A_data :/ WH_data
        div = editvalue(div, 0, 1)
        
        // Reshape div to allow matrix multiplication
        log_div = colshape(log(div), 1)

        res = A_data * log_div

        divergence = res + sum_WH - sum(A_data)

    }
    else if (beta == 2) {
        // 2 = Frobenius (or Euclidean) norm
        delta = A - W*H
        divergence = sqrt(sum(delta :^ 2))
    }
    else {
        _error("Invalid value for beta supplied.")
    }

    return(divergence)
}

void multiplicativeUpdating(real matrix A, 
                            real scalar k, 
                            real scalar iter, 
                            real matrix W, 
                            real matrix H, 
                            real scalar e, 
                            real matrix norms,
                            real scalar beta,
                            real scalar stop)
{
    for (i = 1; i <= iter; i++) {
        
        // Update H
        real matrix W_TA, W_TWH
        W_TA = W' * A
        W_TWH = W' * W * H + J(rows(H), cols(H), e)
        H = H :* W_TA :/ W_TWH

        // Update W
        real matrix AH_T, WHH_T
        AH_T = A * H'
        WHH_T = W * H * H' + J(rows(W), cols(W), e)
        W = W :* AH_T :/ WHH_T
        
        // Calculate norm of difference matrix
        normResult = betaDivergence(A, W, H, beta)

        // Update norms matrix with result of current iteration
        if (i == 1) {
            norms[1, 1] = normResult
        }
        else {
            norm = J(1, 1, normResult)
            norms = norms \ norm
        }

        // Print result of iteration to the screen
        betaMethodString = ""
        if (beta == 0) {
            betaMethodString = "Itakura-Saito Divergence"
        }
        else if (beta == 1) {
            betaMethodString =  "Generalized Kullback-Leibler Divergence"
        }
        else if (beta == 2) {
            betaMethodString = "Frobenius (Euclidean) Norm"
        }

        if (i == 1 | mod(i, 10) == 0 | i == iter) {
            printf("Iteration " + strofreal(i) + " of " + strofreal(iter) + ":\t\tLoss - " + betaMethodString + ":  %9.2f\n", normResult)
        }
        
        // Implement stopping rule if one is set (i.e. stop > 0)
        // Checks every 10 iterations whether MU should contine or if it should stop.
        if (stop != 0 & mod(i, 10) == 0){

            // if ((previous error - current error) / error at initiation) < stop tolerance) then stop
            stopMeasure = (norms[i - 1, 1] - norms[i, 1]) / norms[1, 1]
            if (stopMeasure < stop)
            {
                printf("\nStopping at iteration " + strofreal(i) + "...\n")
                printf("Criteria for early stoping have been met. Error reduced to: " + strofreal(stopMeasure) + ", which < the stopping threshold (" + strofreal(stop) +").\n\n")
                break
            }
        }
    }
}

void coordinateDescentUpdating(real matrix A, 
                               real scalar k, 
                               real scalar iter, 
                               real matrix W, 
                               real matrix H, 
                               real scalar e, 
                               real matrix norms,
                               real scalar beta,
                               real scalar stop)
{
    for (i = 1; i <= iter; i++) {
        

        //http://www.cs.cornell.edu/~bindel/class/sjtu-summer18/lec/2018-06-21.pdf
        //https://gist.github.com/mblondel/09648344984565f9477a

        // line 420: https://github.com/scikit-learn/scikit-learn/blob/36958fb240fbe435673a9e3c52e769f01f36bec0/sklearn/decomposition/_nmf.py#L420


    }
}



end