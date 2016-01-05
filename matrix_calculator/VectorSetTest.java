package matrix_calculator;

import static org.junit.Assert.*;

import org.junit.Test;

public class VectorSetTest {
    @Test
    public void test() {
        Vector vector = new Vector(1, 2, 3, 4);
        assertEquals(Math.sqrt(30), vector.magnitude(), 0.1); // magnitude works
        assertEquals(4, vector.numRows()); // numrows works
        double[] normalized = vector.normalize();
        assertEquals(1.0 / Math.sqrt(30), normalized[0], 0);
        assertEquals(2.0 / Math.sqrt(30), normalized[1], 0);
        assertEquals(3.0 / Math.sqrt(30), normalized[2], 0);
        assertEquals(4.0 / Math.sqrt(30), normalized[3], 0); // normalization works
        VectorSet set = new VectorSet(vector);
        assertEquals(set.size(), 1);
        set.gramSchmidt(0, null);
        assertEquals(vector, set.get(0)); // gram-schmidt for single vector
                                          // works
        Vector vector1 = new Vector(1, 2, 3, 0);
        Vector vector2 = new Vector(1, 2, 0, 0);
        Vector vector3 = new Vector(1, 0, 0, 1);
        set = new VectorSet(vector1, vector2, vector3);
        set = set.gramSchmidt(0, null);
        set.print();
        assertEquals(set.get(0).values()[2], 3, 0);
        assertEquals(set.get(1).values()[0], 9.0 / 14.0, 0.1);
        assertEquals(set.get(1).values()[1], 9.0 / 7.0, 0.1);
        assertEquals(set.get(1).values()[2], -15.0 / 14.0, 0.1);
        assertEquals(set.get(1).values()[3], 0, 0.1);
        assertEquals(set.get(2).values()[0], 4.0 / 5.0, 0.0);
        assertEquals(set.get(2).values()[1], -2.0 / 5.0, 0.0);
        assertEquals(set.get(2).values()[2], 0, 0.0000000001);
        assertEquals(set.get(2).values()[3], 1, 0.0000000001);
    }

    @Test
    public void testBases() throws MatrixException {
        Vector vector1 = new Vector(1, 0);
        vector1.matricize().print();
        Vector vector2 = new Vector(1, 2);
        vector2.matricize().print();
        Vector coordinateVector = new Vector(-2, 3);
        VectorSet set = new VectorSet(vector1, vector2);
        set.matricize().print();
        Vector coordinates = set.coordinates(coordinateVector);
        assertEquals(coordinates.values()[0], 1, 0);
        assertEquals(coordinates.values()[1], 6, 0); // Finding coordinates
        // works
        vector1 = new Vector(-9, 1);
        vector2 = new Vector(-5, -1);
        Vector vector3 = new Vector(1, -4);
        Vector vector4 = new Vector(3, -5);
        VectorSet bSet = new VectorSet(vector1, vector2);
        VectorSet cSet = new VectorSet(vector3, vector4);
        Matrix bMatrix = bSet.matricize();
        Matrix cMatrix = cSet.matricize();
        bMatrix.changeOfBasis(cMatrix).print();
        vector1 = new Vector(6, 4);
        vector2 = new Vector(4, 1);
        vector3 = new Vector(-6, 1);
        set = new VectorSet(vector2, vector3);
        set.print();
        SquareMatrix matrix1 = (SquareMatrix) set.matricize();
        Matrix matrix2 = matrix1.getInverse();
        System.out.println("lol");
        System.out.println("ok");
        Operations.matrixMult(matrix1, matrix2).print();
    }

    
    @Test
    public void testEquality() throws MatrixException {
        Vector vector1 = new Vector(1, 2);
        Vector vector2 = new Vector(4, 3);
        Vector vector3 = new Vector(5, 6);
        VectorSet set1 = new VectorSet(vector1, vector2, vector3);
        Vector vector4 = new Vector(4, 3);
        Vector vector5 = new Vector(5, 6);
        Vector vector6 = new Vector(1, 2);
        VectorSet set2 = new VectorSet(vector4, vector5, vector6);
        assertTrue(vector1.equals(vector4));
        assertTrue(set1.contains(vector4));
        assertTrue(set1.contains(vector5));
        assertTrue(set1.contains(vector6));
        assertTrue(set1.equals(set2));
    }
    
    @Test
    public void testColSpace() throws MatrixException {
        double[][] contents = { { -2.0, -5.0, 8.0, 0, -17 }, { 1.0, 3.0, -5.0, 1.0, 5.0 },
                { 3.0, 11.0, -19.0, 7.0, 1.0 }, { 1.0, 7.0, -13.0, 5.0, -3.0 } };
        Matrix matrix = new Matrix(4, 5, contents);
        VectorSet basis = matrix.columnSpace();
        basis.print();
    }
}
