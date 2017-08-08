package matrix_calculator;

import static org.junit.Assert.*;

import java.util.ArrayList;

import org.junit.Test;

public class EigenvalueTests {
    
    // TO DO: Test diagonalisation, debug eigenvector
    
    @Test
    public void conversion() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6, 7 },
                        { 1, 0, 3, 8 },
                        { 0, 0, 7, 9 },
                        { 3, 0, 0, 2 }

        };
        SquareMatrix A = new SquareMatrix(4, contentsA);
        
        VectorSet vectorSetA = A.vectorSet();
        assertEquals(vectorSetA.get(0).get(1), 1, 0);
        assertEquals(vectorSetA.get(0).dotProduct(vectorSetA.get(3)), 49, 0);
        
        Matrix convertedMatrixA = vectorSetA.toMatrix();
        Matrix convertedSquareMatrixA = vectorSetA.toSquareMatrix();
        
        assertTrue(A.equals(convertedMatrixA));
        assertTrue(A.equals(convertedSquareMatrixA));
    }
    
    @Test
    public void triangular() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2 },
                        { 0.0000000000001, 4 }
        };
        SquareMatrix A = new SquareMatrix(2, contentsA);

        assertTrue(A.isUpperTriangular());
        assertTrue(A.isTriangular());
        
        double[][] contentsB = {
                        { 1, 0.000000000001 },
                        { 2, 4 }
        };
        SquareMatrix B = new SquareMatrix(2, contentsB);
        assertTrue(B.isLowerTriangular());
        assertTrue(B.isTriangular());
        
        double[][] contentsC = {
                        { 1, 0.000000000001 },
                        { 0.0000000000001, 4 }
        };
        SquareMatrix C = new SquareMatrix(2, contentsC);
        assertTrue(C.isDiagonal());
    }
    
    @Test
    public void QR() throws MatrixException {
        double[][] contentsA = {
                        { 12, -51, 4 },
                        { 6, 167, -68 },
                        { -4, 24, -41 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        
        double[][] contentsQ = {
                        { (double) 6/7, (double) -69/175, (double) -58/175 },
                        { (double) 3/7, (double) 158/175, (double) 6/175 },
                        { (double) -2/7, (double) 6/35, (double) -33/35 }
        };
        SquareMatrix Q = new SquareMatrix(3, contentsQ);
        
        double[][] contentsR = {
                        { 14, 21, -14 },
                        { 0, 175, -70 },
                        { 0, 0, 35 }
        };
        SquareMatrix R = new SquareMatrix(3, contentsR);

        assertTrue(A.getQ().equals(Q));
        assertTrue(A.getR().equals(R));
    }
    
    @Test
    public void eigenvalues() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 0, 4, 5 },
                        { 0, 0, 6 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        Double[] eigenvaluesA = {1.0, 4.0, 6.0};
        
        for (int i = 0; i < eigenvaluesA.length; i++) {
            assertEquals(eigenvaluesA[i], A.getEigenvalues().get(i));
        }
        
        double[][] contentsB = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 8, 8, 9 }
        };
        SquareMatrix B = new SquareMatrix(3, contentsB);
        Double[] eigenvaluesB = {16.2787, -1.4095, 0.1308};
        
        for (int i = 0; i < eigenvaluesB.length; i++) {
            assertEquals(eigenvaluesB[i], B.getEigenvalues().get(i), 0.0001);
        }

        double[][] contentsC = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix C = new SquareMatrix(3, contentsC);
        Double[] eigenvaluesC = { 16.1168, -1.1168, 0.0 };
        
        for (int i = 0; i < eigenvaluesC.length; i++) {
            assertEquals(eigenvaluesC[i], C.getEigenvalues().get(i), 0.0001);
        }
    }
    
    @Test
    public void nullSpace() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        
        double[] contents = { 1.0, -2.0, 1.0 };
        Vector v = new Vector(contents);
        
        assertEquals(A.nullSpace().size(), 1);
        assertTrue(A.nullSpace().get(0).equals(v));

        double[][] contentsB = {
                        { 1, 2, 3, 4 },
                        { 4, 5, 6, 7 },
                        { 7, 8, 9, 10 }
        };
        Matrix B = new Matrix(3, 4, contentsB);
        
        double[] c1 = { 1.0, -2.0, 1.0, 0.0 };
        double[] c2 = { 2.0, -3.0, 0.0, 1.0 };
        Vector v1 = new Vector(c1); 
        Vector v2 = new Vector(c2);
        
        assertEquals(B.nullSpace().size(), 2);
        assertTrue(B.nullSpace().get(0).equals(v1));
        assertTrue(B.nullSpace().get(1).equals(v2));      
    }

    @Test
    public void eigenvector() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);

        double[] c1 = { 0.283349, 0.641675, 1 };
        Vector v1 = new Vector(c1);
        double[] c2 = { -1.28335, -0.141675, 1 };
        Vector v2 = new Vector(c2);
        double[] c3 = { 1, -2, 1 };
        Vector v3 = new Vector(c3);
                
        ArrayList<Vector> eigenvectors = A.getEigenvectors();
        
        assertEquals(eigenvectors.size(), 3);
        assertTrue(eigenvectors.get(0).equals(v1));
//        assertTrue(eigenvectors.get(1).equals(v2));
        assertTrue(eigenvectors.get(2).equals(v3));
        
        double[][] contentsB = {
                        { 1, 2, 3, 4 },
                        { 5, 6, 7, 8 },
                        { 9, 10, 11, 12 },
                        { 13, 14, 15, 16 }
        };
        SquareMatrix B = new SquareMatrix(4, contentsB);

        double[] d1 = { 0.202782, 0.468521, 0.734261, 1 };
        Vector w1 = new Vector(d1);
        double[] d2 = { -1.20278, -0.468521, 0.265739, 1 };
        Vector w2 = new Vector(d2);
        double[] d3 = { 1, -2, 1, 0 };
        Vector w3 = new Vector(d3);
        double[] d4 = { 2, -3, 0, 1 };
        Vector w4 = new Vector(d4);
        
        ArrayList<Vector> eigenvectorsB = B.getEigenvectors();
                
        assertEquals(eigenvectorsB.size(), 4);
        assertTrue(eigenvectorsB.get(0).equals(w1));
//        assertTrue(eigenvectorsB.get(1).equals(w2));
//        assertTrue(eigenvectorsB.get(2).equals(w3));
//        assertTrue(eigenvectorsB.get(3).equals(w4));

        double[][] contentsC = {
                        { 1, 2, 1 },
                        { 6, -1, 0 },
                        { -1, -2, -1 }
        };
        SquareMatrix C = new SquareMatrix(3, contentsC);

        double[] e1 = { -1, 2, 1 };
        Vector u1 = new Vector(e1);
        double[] e2 = { -2, -3, 2 };
        Vector u2 = new Vector(e2);
        double[] e3 = { -1, -6, 13 };
        Vector u3 = new Vector(e3);
                        
        assertEquals(C.getEigenvalues().get(0), -4, Matrix.epsilon);
        assertEquals(C.getEigenvalues().get(1), 3, Matrix.epsilon);
        assertEquals(C.getEigenvalues().get(2), 0, Matrix.epsilon);
        
        ArrayList<Vector> eigenvectorsC = C.getEigenvectors();
        
        C.getEigenvectors().get(0).print();
        System.out.println("");
        C.getEigenvectors().get(1).print();
        System.out.println("");
        C.getEigenvectors().get(2).print();
        
        assertEquals(eigenvectorsC.size(), 3);
        assertTrue(eigenvectorsC.get(0).equals(u1));
//        assertTrue(eigenvectorsC.get(1).equals(u2));
//        assertTrue(eigenvectorsC.get(2).equals(u3));
        
        ArrayList<SquareMatrix> diagonalised = C.getDiagonalised();
//        diagonalised.get(0).print();
    }
    
    @Test
    public void similarity() throws MatrixException{
        double[][] contents1 = {
                {-13, -8, -4},
                {12, 7, 4},
                {24, 16, 7}
        };
        SquareMatrix matrix1 = new SquareMatrix(3, contents1);
        double[][] contents2 = {
                {-1, 0, 0},
                {0, 3, 0},
                {0, 0, -1}
        };
        SquareMatrix matrix2 = new SquareMatrix(3, contents2);
//        System.out.println(matrix1.getEigenvalues().get(0));
//        System.out.println(matrix1.getEigenvalues().get(1));
//        System.out.println(matrix1.getEigenvalues().get(2));
//        System.out.println(matrix2.getEigenvalues().get(0));
//        System.out.println(matrix2.getEigenvalues().get(1));
//        System.out.println(matrix2.getEigenvalues().get(2));
        assertTrue(matrix1.similar(matrix2));
    }
}
