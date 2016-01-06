package matrix_calculator;

import static org.junit.Assert.*;

import org.junit.Test;

public class SolutionTest {

    @Test
    public void linearsystem() throws MatrixException {
        double[][] contentsA = {
                        { 1, 2, 3 },
                        { 4, 5, 6 },
                        { 7, 8, 9 }
        };
        SquareMatrix A = new SquareMatrix(3, contentsA);
        
        double[] c1 = {12, 16, 9};
        Vector b1 = new Vector(c1);
        
        Vector x1 = A.solve(b1); // Should have no solution
        assertEquals(x1, null);
        
        double[][] contentsB = {
                        { 1, 0, 1, 0 },
                        { 0, 2, 2, 2 },
                        { 4, -2, 2, -2 }
        };
        Matrix B = new Matrix(3, 4, contentsB);
        
        double c2[] = {2, -10, 18};
        Vector b2 = new Vector(c2);
        
        Vector x2 = B.solve(b2);
        
        double[] r2 = {2.0, -5.0, 0.0, 0.0};
        Vector result2 = new Vector(r2);
        
        assertTrue(x2.equals(result2));
        
        B.generalSolution(b2);
    }

    @Test
    public void axB() throws MatrixException {
        double[][] contents = {
                        { 3, 5, -4 },
                        { -3, -2, 4 },
                        { 6, 1, -8 }
        };
        Matrix matrix = new Matrix(3, 3, contents);
        Vector vector = new Vector(7, -1, -4);
        Vector solution = new Vector(-1, 2, 0);
        assertTrue(matrix.solve(vector).equalScale(solution));

        double[][] contents2 = {
                        { 1, 2, -1 },
                        { 0, -5, 3 },

        };
        Matrix matrix2 = new Matrix(2, 3, contents2);
        Vector vector2 = new Vector(3, 6);
        Vector solution2 = new Vector(4, 3, 7);

        matrix2.generalSolution(vector2);

        double[][] contents3 = {
                        { 1, 3, 0, 2 },
                        { 0, 0, 1, 4 },
                        { 1, 3, 1, 6 }
        };
        Matrix matrix3 = new Matrix(3, 4, contents3);
        Vector vector3 = new Vector(1, 6, 7);
        Vector solution3 = new Vector(1, 0, 6, 0);
        Vector answer = matrix3.solve(vector3);
        assertTrue(answer.equalScale(solution3));
    }
    
}
