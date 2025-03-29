static void main(String[] args) {
    // Example matrix 0 - 2x2 matrix
//        def matrix = [
//                [7, 3],
//                [3, -1]
//        ]

    // Example matrix 1 - 2x2 matrix
    def matrix = [
            [4, -5],
            [2, -3]
    ]

    //Example matrix 2 - 2x2 matrix with complex eigenvalues

//        def matrix =[
//                [0, -1],
//                [1,0]
//        ]

    //Example matrix 3 - 3x3 matrix

//        def matrix = [
//                [2, 2, -2],
//                [1, 3, -1],
//                [-1, 1, 1],
//        ]

    // Example matrix 4: 3x3 triangular matrix
//        def matrix = [
//                [1, 2,4],
//                [0, 4, 7],
//                [0, 0, 6],
//        ]


    println "Original Matrix:"
    printMatrix(matrix)

    // Step 1: Find eigenvalues
    def eigenvalues = findEigenvalues(matrix)
    println "\nEigenvalues:"
    eigenvalues.eachWithIndex { val, i ->
        println "λ${i + 1} = ${val}"
    }

    // Step 2: Find eigenvectors and eigenspaces
    println "\nEigenvectors and Eigenspaces:"
    eigenvalues.eachWithIndex { eigenvalue, i ->
        def eigenspace = findEigenspace(matrix, eigenvalue)
        println "\nFor eigenvalue λ${i + 1} = ${eigenvalue}:"
        println "Eigenvectors: ${eigenspace.basis}"
        println "Eigenspace: span${eigenspace.basis}"
    }
}

// Finding eigenvalues
static def findEigenvalues(List<List<Double>> matrix) {
    if (matrix.size() == 2 && matrix[0].size() == 2) {
        return getSolutionOfQuadratic(matrix)
    } else if (matrix.size() == 3 && matrix[0].size() == 3) {
        return getSolutionOfCubic(matrix)
    }
    throw new IllegalArgumentException("Matrix must be quadratic or cubic (2x2 or 3x3)")
}

//  [ [a, b], [c, d] ]
static List<Double> getSolutionOfQuadratic(List<List<Double>> m) {
    def a = 1.0
    def b = -(m[0][0] + m[1][1])
    def c = m[0][0] * m[1][1] - m[0][1] * m[1][0]

    def discriminant = b * b - 4 * a * c

    // Check for complex eigenvalues
    if (discriminant < 0) {
        throw new ArithmeticException("Matrix must be real eigenvalues, complex eigenvalues not supported in this implementation.")
    }

    def eigenvalueOne = (-b + Math.sqrt(discriminant)) / (2.0 * a)
    def eigenvalueTwo = (-b - Math.sqrt(discriminant)) / (2.0 * a)
    return [eigenvalueOne, eigenvalueTwo]
}

// Solves eigenvalues for 3x3 matrix
static List<Double> getSolutionOfCubic(List<List<Double>> m) {
    // First check if matrix is triangular
    if (isTriangular(m)) {
        // For triangular matrices, return diagonal elements
        println("Matrix is triangular, eigenvalues are diagonal elements.")
        return (0..<m.size()).collect { m[it][it] }
    }
    def a = m[0][0], b = m[0][1], c = m[0][2]
    def d = m[1][0], e = m[1][1], f = m[1][2]
    def g = m[2][0], h = m[2][1], i = m[2][2]

    // Coefficients for characteristic equation: λ³ + Aλ² + Bλ + C = 0
    def A = -(a + e + i)
    def B = a * e + a * i + e * i - b * d - c * g - f * h
    def C = -(a * e * i + b * f * g + c * d * h - a * f * h - b * d * i - c * e * g)

    // Find roots numerically
    return findRealRoots(1.0, A, B, C)
}

static isTriangular(List<List<Double>> m) {
    def n = m.size()
    boolean isUpper = true
    boolean isLower = true

    // upper triangular
    (1..<n).each { i ->
        (0..<i).each { j ->
            if (m[i][j] != 0) {
                isUpper = false
                return // break out of loop
            }
        }
        if (!isUpper) return
    }

    // lower triangular (zeros above diagonal)
    (0..<n - 1).each { i ->
        (i + 1..<n).each { j ->
            if (m[i][j] != 0) {
                isLower = false
                return // break out of loop
            }
        }
        if (!isLower) return
    }

    return isUpper || isLower
}

// Find real roots of cubic equation
static List<Double> findRealRoots(double a, double b, double c, double d) {
    def roots = []
    double step = 0.01
    double range = 10.0
    double tolerance = 1e-6

    // Convert to int steps to avoid BigDecimal issue
    int numSteps = (range * 2 / step) as int

    (0..numSteps).each { int i ->
        double x = -range + (i * step)
        double y = a * x * x * x + b * x * x + c * x + d

        if (Math.abs(y) < tolerance) {
            if (roots.every { Math.abs(it - x) > step * 2 }) {
                roots << x
            }
        }
    }

    if (roots.isEmpty()) {
        throw new ArithmeticException("No real roots found - complex eigenvalues exist")
    }
    return roots
}

// Finds eigenspace
static Map findEigenspace(List<List<Double>> matrix, Double eigenvalue) {
    // Create A - λI matrix
    def n = matrix.size()
    def aMinus = (0..<n).collect { i ->
        (0..<n).collect { j ->
            matrix[i][j] - (i == j ? eigenvalue : 0)
        }
    }

    // Find null space (A-λI)x = 0
    def basis = findNullSpace(aMinus)

    return [eigenvalue: eigenvalue, basis: basis]
}

// Finds null space of a matrix (simplified implementation)
static List<List<Double>> findNullSpace(List<List<Double>> matrix) {
    // This is a simplified implementation that works for 2x2 and simple cases
    def n = matrix.size()

    if (n == 2) {
        def a = matrix[0][0], b = matrix[0][1]
        def c = matrix[1][0], d = matrix[1][1]

        // Handle zero matrix case
        if ([a, b, c, d].every { it == 0 }) {
            return [[1.0, 0.0], [0.0, 1.0]]
        }

        // If rows are linearly dependent
        if (Math.abs(a * d - b * c) < 1e-10) {
            if (b != 0) {
                def gcd = computeGCD(Math.abs(b), Math.abs(a))
                return [[b / gcd, -a / gcd]]
            } else if (d != 0) {
                return [[d, -c]]
            } else {
                return [[0.0, 1.0]]
            }
        }
    }

    // For 3x3 matrices
    if (n == 3) {
        // Real implementation would need Gaussian elimination
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    }

    //  default vector
    return [[1] * n] // fallback - returns vector of 1s
}

static int computeGCD(double a, double b) {
    return b == 0 ? a : computeGCD(b, a % b)
}

static void printMatrix(matrix) {
    matrix.each { row ->
        println row.join("\t")
    }
}

