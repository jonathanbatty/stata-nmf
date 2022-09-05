capture prog drop nmf
*! version 0.01 05 Sep 2022
program define nmf, rclass
    version 17
 
    syntax varlist(numeric), k(integer) iter(integer) [initial(string) beta(numlist integer max = 1) nograph noframes]
    
    // The aim of NMF is to find W (u x k) and H (k x v) such that A ~ WH, where all 3 matrices contain only non-negative values
    // A good appoximation may be achieved with k << rank(A)

    // Initialisation options, specified using the init() parameter, include:
    //    random [*default] - random initialisation of matrix in range [0, 1] 
    //    nndsvd - Nonnegative Double Singular Value Decomposition [Boutsidis2007] - suited to sparse factors
    //    nndsvda - NNSVD with zero elements replaced with the input matrix average (not recommended for sparse data)
    //    nndsvdar - NNSVD with zero elements replaced with a random value (nor recommended for sparse data)

    // Beta divergence options, specified using the beta() parameter, include:
    //    0 - Itakura-Saito divergence
    //    1 - Generalized Kullback-Leibler divergence
    //    2 [*default] - Frobenius (Euclidean) norm

    // nograph suppresses production of the graph of beta divergence at each iteration
    // noframes does not save frames containing output matrices: W, H and norms 

    // Uses multiplicative update method of Lee and Seung in 1999: 
    // https://proceedings.neurips.cc/paper/2000/file/f9d1152547c0bde01830b7e8bd60024c-Paper.pdf

    // To do: 
    // + Add a stopping condition, e.g. stoptolerance with default 1.0e-4 (like in scikit learn: https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.NMF.html and https://github.com/scikit-learn/scikit-learn/blob/36958fb24/sklearn/decomposition/_nmf.py#L1158)
    /// -> e.g. for all iterations beyond first, compare norms with previous value (i - 1)... if two successive norms are within a certain threshold; stop
    // Consider scaling random values -? by sqrt(average of input matrix / length of input matrix) as per skl?
    // Test and benchmark
    // 

    // If a value for initialisation method is not pased, default option is is random initialisation 
    if "`initial'" == "" local initial "random"

    // If a value specifying the form of beta divergence normalisation is not pased, default option is is Frobenius normalisation (2) 
    if "`beta'" == "" local beta = 2

    // Run NMF
    mata: nmf("`varlist'", `k', `iter', "`initial'", `beta')
    
    // Store using stata matrices
    matrix W = r(W)
    matrix H = r(H)
    matrix norms = r(norms)

    // Creates frames containing output matrices: W, H and norms
    if "`frames'" != "noframes" {

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
        frame norms: graph twoway line norm iteration, yscale(range(0 .)) ylabel(#5)
    }

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
         real scalar beta)
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
    if (k < 2 | k >= cols(A)) {
        _error("The rank of the output matrix must be at least 2, but less than the rank of the input matrix, A.")
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
    norms = J(iter, 1, .)

    // Perform multiplicative updating given number of iterations
    multiplicativeUpdating(A, k, iter, W, H, e, norms, beta)
        
    // Return results object to stata
    st_matrix("r(W)", W)
    st_matrix("r(H)", H)
    st_matrix("r(norms)", norms)

    return(nmf)
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
        divergence = norm(delta, 0)
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
                            real scalar beta)
{
    for (i = 1; i <= iter; i++) {
        
        printf("Iteration " + strofreal(i) + " of " + strofreal(iter) + "\n")

        // Update H
        W_TA = W' * A
        W_TWH = W' * W * H + J(rows(H), cols(H), e)
        H = H :* W_TA :/ W_TWH
        
        // Update W
        AH_T = A * H'
        WHH_T = W * H * H' + J(rows(W), cols(W), e)
        W = W :* AH_T :/ WHH_T
        
        // Calculate norm of difference matrix
        normResult = betaDivergence(A, W, H, beta)

        // Update norms matrix with result of current iteration
        norms[i, 1] = normResult
    }
}

end