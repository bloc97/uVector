/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vectors;

/**
 *
 * @author bowen
 */
public class Vector4 implements Vector<Vector4> {

    private final double[] content;
    
    public Vector4() {
        content = new double[4];
    }

    @Override
    public double[] content() {
        return content;
    }

    @Override
    public Vector4 set(int i, double d) {
        content[i] = d;
        return this;
    }

    @Override
    public Vector4 shell() {
        return new Vector4();
    }
    
}
