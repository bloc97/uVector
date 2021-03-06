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
public class VectorN implements Vector<VectorN> {

    private final double[] content;
    
    public VectorN(int n) {
        content = new double[n];
    }

    @Override
    public double[] content() {
        return content;
    }

    @Override
    public VectorN set(int i, double d) {
        content[i] = d;
        return this;
    }

    @Override
    public VectorN shell() {
        return new VectorN(size());
    }
    
}
