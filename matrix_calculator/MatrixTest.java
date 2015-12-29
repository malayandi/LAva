package matrix_calculator;

import static org.junit.Assert.*;
import org.junit.Test;

public class MatrixTest {

    public void basicFunctionality() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 6, 3, 2 }
        };
        Matrix A = new Matrix(2, 3, contents);

        assertEquals(A.getHeight(), 2);
        assertEquals(A.getWidth(), 3);
        assertEquals(A.get(1, 2), 4, 0);

        Matrix B = new Matrix(2, 3, new double[2][3]);
        B.set(1, 1, 5);
        B.set(1, 2, 4);
        B.set(1, 3, 6);
        B.set(2, 1, 6);
        B.set(2, 2, 3);
        B.set(2, 3, 2);
        assertTrue(A.equals(B));

        A.set(2, 3, 15.9);
        assertEquals(A.get(2, 3), 15.9, 0);
        A.print();        
    }

    public void elementaryRowOperations() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 6, 3, 2 }
        };
        Matrix A = new Matrix(2, 3, contents);

        System.out.println("");
        A.scalarMultRow(1, 3);
        A.print();

        System.out.println("");
        A.scalarMult(2);
        A.print();

        System.out.println("");
        A.switchRow(1, 2);
        A.print();

        System.out.println("");
        A.add(1, 2, 2);
        A.print();
    }

    public void transpose() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contentsA);
        double[][] contentsB = {
                        { 5, 1, 0, 3 },
                        { 4, 0, 0, 0 },
                        { 6, 3, 7, 0 },

        };
        Matrix B = new Matrix(3, 4, contentsB);
        assertTrue(A.getTranspose().equals(B));
    }


    public void count() throws MatrixException {
        double[][] contents = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contents);

        assertEquals(A.count(1), 3);
        assertEquals(A.count(2), 1);
        assertEquals(A.count(3), 3);

        assertEquals(A.count(1, 3), 1);
        assertEquals(A.count(2, 2), 0);
        assertEquals(A.count(2, 3), 0);
    }

    public void basicOperations() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contentsA);
        double[][] contentsB = {
                        { 4, 4, 7 },
                        { 3, 8, 3 },
                        { 7, 12, 7 },
                        { 3, 5, 9 }

        };
        Matrix B = new Matrix(4, 3, contentsB);

        System.out.println("");
        Matrix C = Operations.add(A, B);
        C.print();

        System.out.println("");
        Matrix D = Operations.scalarMult(A, 2.5);
        D.print();

        System.out.println("");
        Matrix E = Operations.subtract(A, B);
        E.print();
    }

    public void multiplication() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contentsA);
        double[][] contentsB = {
                        { 1, 4, 3, 6 },
                        { 8, 5, 2, 2 },
                        { 1, 6, 9, 4 },

        };
        Matrix B = new Matrix(3, 4, contentsB);
        double[][] contentsAB = {
                        { 43, 76, 77, 62 },
                        { 4, 22, 30, 18 },
                        { 7, 42, 63, 28 },
                        { 3, 12, 9, 18 }

        };
        Matrix AB = new Matrix(4, 4, contentsAB);
        Matrix ABTest = Operations.matrixMult(A, B);
        assertTrue(AB.equals(ABTest));
        
        double[][] contentsA2 = {
                        { 2, 5 },
                        { 3, 7 }
        };
        SquareMatrix A2 = new SquareMatrix(2, contentsA2);
        double[][] contentsA2Inverse = {
                        { -7, 5 },
                        { 3, -2 }
        };
        SquareMatrix A2Inverse = new SquareMatrix(2, contentsA2Inverse);
        Matrix ITest = Operations.matrixMult(A2, A2Inverse);
        ITest.print();
    }

    public void rowReduction() throws MatrixException {
        double[][] contentsA = {
                        { 5, 4, 6 },
                        { 1, 0, 3 },
                        { 0, 0, 7 },
                        { 3, 0, 0 }

        };
        Matrix A = new Matrix(4, 3, contentsA);

        double[][] contentsRRA = {
                        { 1, 0.8, 1.2 },
                        { 0, 1, -2.25 },
                        { 0, 0, 1 },
                        { 0, 0, 0 }

        };
        Matrix RRA = new Matrix(4, 3, contentsRRA);

        System.out.println("");
        A.getRowRed().print();

        System.out.println("");
        A.getRowRedEF().print();

        assertTrue(A.isLinInd());
        assertTrue(A.isInjective());
        assertFalse(A.isSurjective());
        assertEquals(A.getRank(), 3);
        assertEquals(A.getNullity(), 0);


        double[][] contentsB = {
                        { 0, 4, 6, 7 },
                        { 0, 2, 3, 12 },
                        { 0, 9, 7, 0 },
        };
        Matrix B = new Matrix(3, 4, contentsB);

        double[][] contentsRRB = {
                        { 0, 1, 1.5, 1.75 },
                        { 0, 0, 1, 63/26 },
                        { 0, 0, 0, 1 },
        };
        Matrix RRB = new Matrix(3, 4, contentsRRB);

        System.out.println("");
        B.getRowRed().print();

        System.out.println("");
        B.getRowRedEF().print();

        assertFalse(B.isLinInd());
        assertFalse(B.isInjective());
        assertTrue(B.isSurjective());
        assertEquals(B.getRank(), 3);
        assertEquals(B.getNullity(), 1);
    }

    @Test
    public void basicSquareMatrix() throws MatrixException {
        double[][] contentsA = {
                        { 1, 5 },
                        { 2, 7 }
        };
        SquareMatrix A = new SquareMatrix(2, contentsA);
        double[][] contentsAT = {
                        { 1, 2 },
                        { 5, 7 }
        };
        SquareMatrix AT = new SquareMatrix(2, contentsAT);
        double[][] contentsAInverse = {
                        { (double) -7/3, (double) 5/3 },
                        { (double) 2/3, (double) -1/3 }
        };
        SquareMatrix AInverse = new SquareMatrix(2, contentsAInverse);
        
        assertTrue(A.getTranspose().equals(AT));
        assertEquals(A.getTrace(), 8 , 0);
        assertEquals(A.getDet(), -3, 0);
        assertEquals(A.getDet(), 1/A.getInverse().getDet(), 0); 
        assertTrue(A.getInverse().equals(AInverse));
    }
}



