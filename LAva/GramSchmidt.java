package matrix_calculator;

import static org.junit.Assert.*;

import org.junit.Test;

public class GramSchmidt {

    @Test
    public void test() {
        Vector v1 = new Vector(1,1,1);
        Vector v2 = new Vector(-1,0,1);
        Vector v3 = new Vector(-1,1,0);
        VectorSet v = new VectorSet(v1,v2,v3);
        v = gramSchmidt(v);
    }

}
