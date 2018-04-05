/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package vectors;

import java.util.function.BiFunction;

/**
 *
 * @author bowen
 */
public interface DistanceMetric extends BiFunction<Vector, Vector, Double> {
    public static DistanceMetric 
            TAXICAB = getLpSpaceMetric(1),
            EUCLIDIAN = getLpSpaceMetric(2),
            MAX = getLpSpaceMetric(Double.POSITIVE_INFINITY),
            RANDOM = (t, u) -> {
        return Math.random();
    };
    
    public static DistanceMetric getLpSpaceMetric(double p) {
        return (t, u) -> {
            return Vectors.normLp(t.copy().sub(u), p);
        };
    }

    @Override
    public Double apply(Vector target, Vector u);
    
}
