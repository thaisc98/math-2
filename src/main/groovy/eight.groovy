
class MatrixEigensolver {

    static void main(String[] args) {
        // Example matrix - replace with your own
//        def matrix = [
//                [4, -5],
//                [2, -3]
//        ]
        def matrix =[
                [3, 6],
                [2,1]
        ]


        println "Original Matrix:"
        printMatrix(matrix)

        // Step 1: Find eigenvalues
        def eigenvalues = findEigenvalues(matrix)
        println "\nEigenvalues:"
        eigenvalues.eachWithIndex { val, i ->
            println "fλ${i+1} = ${val}"
        }

        // Step 2: Find eigenvectors and eigenspaces
        println "\nEigenvectors and Eigenspaces:"
        eigenvalues.each { λ ->
            def eigenspace = findEigenspace(matrix, λ)
            println "\nFor eigenvalue λ = ${λ}:"
            println "Eigenvectors: ${eigenspace.basis}"
            println "Eigenspace: span${eigenspace.basis}"
        }
    }

    // Finds eigenvalues by solving the characteristic polynomial
    static List<Number> findEigenvalues(List<List<Number>> matrix) {
        if (matrix.size() == 2 && matrix[0].size() == 2) {
            return solveQuadratic(matrix)
        } else if (matrix.size() == 3 && matrix[0].size() == 3) {
            return solveCubic(matrix)
        } else {
            throw new UnsupportedOperationException("Only 2x2 and 3x3 matrices are supported in this implementation")
        }
    }

    // Solves eigenvalues for 2x2 matrix: [ [a, b], [c, d] ]
    static List<Number> solveQuadratic(List<List<Number>> m) {
        def a = 1
        def b = -(m[0][0] + m[1][1])
        def c = m[0][0]*m[1][1] - m[0][1]*m[1][0]

        def discriminant = b*b - 4*a*c

        if (discriminant >= 0) {
            def lambda1 = (-b + Math.sqrt(discriminant)) / (2*a)
            def lambda2 = (-b - Math.sqrt(discriminant)) / (2*a)
            return [lambda1, lambda2]
        } else {
            def realPart = -b / (2*a)
            def imagPart = Math.sqrt(-discriminant) / (2*a)
            return [new ComplexNumber(realPart, imagPart),
                    new ComplexNumber(realPart, -imagPart)]
        }
    }

    // Solves eigenvalues for 3x3 matrix (simplified implementation)
    static List<Number> solveCubic(List<List<Number>> m) {
        // This is a simplified implementation that only finds real roots
        // For a complete solution, we'd need to implement the full cubic formula
        def a = m[0][0], b = m[0][1], c = m[0][2]
        def d = m[1][0], e = m[1][1], f = m[1][2]
        def g = m[2][0], h = m[2][1], i = m[2][2]

        // Coefficients for characteristic equation: λ³ + Aλ² + Bλ + C = 0
        def A = -(a + e + i)
        def B = a*e + a*i + e*i - b*d - c*g - f*h
        def C = -(a*e*i + b*f*g + c*d*h - a*f*h - b*d*i - c*e*g)

        // Find roots numerically (for simplicity)
        return findRealRoots(1, A, B, C)
    }

    // Numerically find real roots of cubic equation (simplified)
    static List<Number> findRealRoots(double a, double b, double c, double d) {
        // This is a very basic root finder - in practice you'd want a more robust method
        def roots = []
        def step = 0.1
        def range = 10.0

        (-range..range).step(step) { x ->
            def y = a*x*x*x + b*x*x + c*x + d
            if (Math.abs(y) < 0.1) { // crude threshold
                if (roots.every { Math.abs(it - x) > 0.5 }) { // avoid duplicates
                    roots << x
                }
            }
        }

        return roots.size() > 0 ? roots : [0] // fallback
    }

    // Finds eigenspace for a given eigenvalue
    static Map findEigenspace(List<List<Number>> matrix, Number λ) {
        // Create A - λI matrix
        def n = matrix.size()
        def aMinus = (0..<n).collect { i ->
            (0..<n).collect { j ->
                matrix[i][j] - (i == j ? λ : 0)
            }
        }

        // Find null space (solutions to (A-λI)x = 0)
        def basis = findNullSpace(aMinus)

        return [eigenvalue: λ, basis: basis]
    }

    // Finds null space of a matrix (simplified implementation)
    static List<List<Number>> findNullSpace(List<List<Number>> matrix) {
        // This is a simplified implementation that works for 2x2 and simple cases
        def n = matrix.size()

        if (n == 2) {
            def a = matrix[0][0], b = matrix[0][1]
            def c = matrix[1][0], d = matrix[1][1]

            // Handle zero matrix case
            if ([a,b,c,d].every { it == 0 }) {
                return [[1, 0], [0, 1]]
            }

            // If one row is multiple of another
            if (a != 0 && c != 0 && b/a == d/c) {
                def ratio = c/a
                if (b != 0) {
                    return [[-b/a, 1]]
                } else {
                    return [[-d/c, 1]]
                }
            }

            // If one column is multiple of another
            if (a != 0 && b != 0 && c/a == d/b) {
                return [[-a/b, 1]]
            }

            // General solution for 2x2
            if (a*d - b*c == 0) { // determinant is zero
                if (b != 0) {
                    return [[-b/a, 1]]
                } else if (d != 0) {
                    return [[-d/c, 1]]
                } else {
                    return [[1, 0]]
                }
            }
        }

        // For 3x3 or more complex cases, we would need proper Gaussian elimination
        // This is a simplified version that returns a default vector
        return [[1] * n] // fallback - returns vector of 1s
    }

    static void printMatrix(matrix) {
        matrix.each { row ->
            println row.collect { it instanceof ComplexNumber ? it.toString() : it }.join("\t")
        }
    }
}

// Helper class for complex numbers
class ComplexNumber {
    double real
    double imaginary

    ComplexNumber(double real, double imaginary = 0) {
        this.real = real
        this.imaginary = imaginary
    }

    Number minus(Number n) {
        new ComplexNumber(real - n, imaginary)
    }

    String toString() {
        if (imaginary == 0) return real.toString()
        if (real == 0) return "${imaginary}i"
        "${real} + ${imaginary}i"
    }

    boolean equals(Object other) {
        if (other instanceof Number) return real == other && imaginary == 0
        if (other instanceof ComplexNumber) {
            return real == other.real && imaginary == other.imaginary
        }
        false
    }
}