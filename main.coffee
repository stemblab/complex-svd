nm = numeric #;

# Helper functions

# Stack arrays horizontally
hstack = (A,B) ->
    r.concat(B[idx]) for r, idx in A

# Get even indexed elements from a row
get_evens = (x) ->
    e for e, c in x when c % 2 is 0

# Matrix size
size = (A) -> [A.length, A[0].length]

# Complex SVD (U*S*V')

csvd = (C) ->

    # Equivalent real matrix
    K = hstack(C.x, nm.mul(C.y,-1)).concat(hstack(C.y, C.x))

    # Numericjs SVD needs more rows than cols
    [m, n] = size K
    shortK = m < n
    KK = if shortK then nm.transpose(K) else K

    # SVD of equivalent real matrix
    svd2 = nm.svd(KK)
    [m, n] = size KK

    # Initialize complex SVD
    Z = -> new nm.T [0], [0]  # Zero vector
    svd1 = U: Z(), S: Z(), V: Z()

    # real SVD -> complex SVD
    
    # real part from top of even indexed cols; imag part from bottom
    setSvd = (svd1, svd2, f, N) ->
        for r in [0..N/2-1]
            svd1[f].x[r] = get_evens(svd2[f][r])
            svd1[f].y[r] = get_evens(svd2[f][r+N/2])        
    
    # U
    setSvd svd1, svd2, "U", m

    # S
    svd1.S.x = nm.diag(get_evens(svd2.S))
    svd1.S.y = svd1.S.x.mul(0) # S is real

    # V
    setSvd svd1, svd2, "V", n

    # Transpose output if originally transposed input at beginning
    svd =
        V: if shortK then svd1.U else svd1.V
        S: nm.getDiag(svd1.S.x) # S is real
        U: if shortK then svd1.V else svd1.U
    
# Complex input
C = new nm.T([[1, 0, 1],[1, 3, 0]],[[0, 0, 0],[1, 0, 0]])       

svd = csvd(C) #;

# Reconstructed input
C2 = svd.U.dot(nm.diag(svd.S)).dot(svd.V.transjugate())

# Export to other blabs
$blab.csvd = csvd
