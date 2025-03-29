static void main(String[] args) {

//    def mat = [
//            [1, 9],
//            [3, 2]
//    ]


      def mat = [
            [ 1, 0, 2],
            [ 3, 0, 0],
            [ 2, 1, 4],
    ]

//    def mat = [
//            [ 1, 0, 2, -1 ],
//            [ 3, 0, 0, 5 ],
//            [ 2, 1, 4, -3 ],
//            [ 1, 0, 5, 0 ]
//    ]

//  def mat = [
//          [0, 0, 0],
//          [0, 0, 0],
//          [0, 0, 0]
//  ]

  println("Determinant of the matrix is:  ${getDeterminant(mat, mat.size())} ")
}

def static getDeterminant(matrix, size) {

  // Validate that the matrix is square
  if (matrix.any { it.size() != size }) {
    throw new IllegalArgumentException("Matrix must be square")
  }

  if (size == 1) {
    return matrix[0][0]
  }

  if (size == 2) {
    return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0]
  }

  def rest = 0;

  for (col in 0..<size) {
    def subMatrix = createSubMatrix(matrix, size, col)
    def sign = col % 2 == 0 ? 1 : -1
    rest += sign * matrix[0][col] * getDeterminant(subMatrix, size - 1)
  }
  return rest
}


def static createSubMatrix(matrix, size, excludeCol) {
  def subMatrix = new int[size - 1][size - 1]
  for (i in 1..<size) {
    int subCol = 0
    for (j in 0..<size) {
      if (j == excludeCol) continue
      subMatrix[i - 1][subCol++] = matrix[i][j]
    }
  }
  return subMatrix
}