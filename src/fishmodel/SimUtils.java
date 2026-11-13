package fishmodel;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;

/**
 * Shared utility methods and constants used across simulation classes.
 * Methods were previously duplicated in every RunSim*.java, FullFarm*.java,
 * TankSimulation, SensorPosSimulationNMBU, and TestEmptyCage file.
 */
public class SimUtils {

    public static final double HYPOXIA_THRESHOLD = 6;

    public static final LinearInterpolator interpol = new LinearInterpolator();

    // -------------------------------------------------------------------------
    // Vertical interpolation
    // -------------------------------------------------------------------------

    /**
     * Set up vertical profile for model grid based on values at set depths.
     * Extrapolates beyond end values. Delegates to the zOffset overload with offset 0.
     *
     * @param res    The array to put interpolated values into
     * @param depths Depths at which values are given (increasing values)
     * @param values Values at given depths
     * @param kmax   Number of vertical layers
     * @param dz     Vertical resolution
     */
    public static void interpolateVertical(double[] res, double[] depths, double[] values, int kmax, double dz) {
        interpolateVertical(res, depths, values, kmax, dz, 0);
    }

    /**
     * Set up vertical profile for model grid based on values at set depths.
     * Extrapolates beyond end values using Apache Commons Math linear interpolation.
     *
     * @param res         The array to put interpolated values into
     * @param depthsOrig  Depths at which values are given (increasing values)
     * @param values      Values at given depths
     * @param kmax        Number of vertical layers
     * @param dz          Vertical resolution
     * @param zOffset     Offset to add to depths before interpolating
     */
    public static void interpolateVertical(double[] res, double[] depthsOrig, double[] values, int kmax, double dz, double zOffset) {
        double[] depths;
        if (Math.abs(zOffset) > 1e-6) {
            depths = new double[depthsOrig.length];
            for (int i = 0; i < depths.length; i++) {
                depths[i] = depthsOrig[i] + zOffset;
            }
        } else {
            depths = depthsOrig;
        }

        double minDepth = depths[0], maxDepth = depths[depths.length - 1],
                topValue = values[0], bottomValue = values[values.length - 1];
        PolynomialSplineFunction interp = interpol.interpolate(depths, values);
        for (int i = 0; i < res.length; i++) {
            double currDepth = ((double) i + 0.5) * dz;
            if (currDepth < minDepth)
                res[i] = topValue;
            else if (currDepth > maxDepth)
                res[i] = bottomValue;
            else
                res[i] = interp.value(currDepth);
        }
    }

    // -------------------------------------------------------------------------
    // Fish affinity field setup
    // -------------------------------------------------------------------------

    /**
     * Sets fish affinity based on a vertical profile, applying zero outside the
     * mask and below fishMaxDepth.
     *
     * @return Sum of o2Affinity values over the cage domain.
     */
    public static double setO2AffinityWithVerticalProfile(int[] cageDims, double dz, double fishMaxDepth,
                                                          boolean[][][] mask, double[] affinityProfile,
                                                          double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;
        for (int i = 0; i < cageDims[0]; i++)
            for (int j = 0; j < cageDims[1]; j++)
                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];
                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;
                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
        return o2AffSum;
    }

    /**
     * Sets fish affinity based on a vertical profile with a radial edge-decrease factor:
     * affinity tapers toward the cage perimeter (outer 50% of radius) using a 4th-power
     * rolloff. Applies zero outside the mask and below fishMaxDepth.
     *
     * @return Sum of o2Affinity values over the cage domain.
     */
    public static double setO2AffinityWithVerticalProfileAndEdgeDecrease(int[] cageDims, double rad, double dxy,
                                                                          double dz, double fishMaxDepth,
                                                                          boolean[][][] mask, double[] affinityProfile,
                                                                          double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;
        for (int i = 0; i < cageDims[0]; i++)
            for (int j = 0; j < cageDims[1]; j++) {
                double xDist = dxy * (double) (i - cageDims[0] / 2);
                double yDist = dxy * (double) (j - cageDims[1] / 2);
                double distFromCenter = Math.sqrt(xDist * xDist + yDist * yDist);

                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        affinity[i][j][k] = affinityProfile[k];
                        o2Affinity[i][j][k] = affinityProfile[k];
                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;
                    }
                    // Decrease affinity towards edge (outer 50% of radius):
                    double edgeFactor = 1.;
                    if (distFromCenter / rad > 0.5) {
                        edgeFactor -= Math.pow((distFromCenter - 0.5 * rad) / (0.5 * rad), 4);
                        affinity[i][j][k] *= edgeFactor;
                        o2Affinity[i][j][k] *= edgeFactor;
                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
            }
        return o2AffSum;
    }

    /**
     * Sets fish affinity based on current O2 concentration: affinity is 1 above
     * avoidanceThresh (4 mg/L) and decreases linearly to 0 at O2 = 0.
     * Applies zero outside the mask and below fishMaxDepth.
     *
     * @return Sum of o2Affinity values over the cage domain.
     */
    public static double setO2AffinityWithAvoidance(int[] cageDims, double dz, double fishMaxDepth,
                                                    boolean[][][] mask, double[][][] o2,
                                                    double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;
        double avoidanceThresh = 4;

        for (int i = 0; i < cageDims[0]; i++)
            for (int j = 0; j < cageDims[1]; j++)
                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        if (o2[i][j][k] > avoidanceThresh) {
                            affinity[i][j][k] = 1;
                            o2Affinity[i][j][k] = 1;
                        } else {
                            affinity[i][j][k] = o2[i][j][k] / avoidanceThresh;
                            o2Affinity[i][j][k] = o2[i][j][k] / avoidanceThresh;
                        }
                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
        return o2AffSum;
    }

    /**
     * Sets fish affinity based on a vertical profile and a preferred direction vector.
     * Fish in the downstream half of the cage (dot product with dirVector > 0) get full
     * affinity; fish in the upstream half get affinity set to zero in o2Affinity.
     * The spatial affinity field is always set to the full profile value.
     *
     * @param dirVector 2-element [dx, dy] preferred direction vector
     * @return Sum of o2Affinity values over the cage domain.
     */
    public static double setO2AffinityWithVerticalProfileAndDirection(int[] cageDims, double dz, double fishMaxDepth,
                                                                      double[] dirVector, boolean[][][] mask,
                                                                      double[] affinityProfile,
                                                                      double[][][] affinity, double[][][] o2Affinity) {
        double o2AffSum = 0;
        double centerX = ((double) cageDims[0]) / 2.0,
                centerY = ((double) cageDims[1]) / 2.0;
        double dirVectorLength = Math.sqrt(dirVector[0] * dirVector[0] + dirVector[1] * dirVector[1]);

        for (int i = 0; i < cageDims[0]; i++)
            for (int j = 0; j < cageDims[1]; j++) {
                double distX = ((double) i) - centerX,
                        distY = ((double) j) - centerY;
                double distTot = Math.sqrt(distX * distX + distY * distY);
                if (distTot > 0) {
                    distX = distX / distTot;
                    distY = distY / distTot;
                }
                double dotProduct = distX * dirVector[0] + distY * dirVector[1];

                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        if (dotProduct / dirVectorLength > 0) {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = affinityProfile[k];
                        } else {
                            affinity[i][j][k] = affinityProfile[k];
                            o2Affinity[i][j][k] = 0;
                        }
                    } else {
                        affinity[i][j][k] = 0;
                        o2Affinity[i][j][k] = 0;
                    }
                    o2AffSum += o2Affinity[i][j][k];
                }
            }
        return o2AffSum;
    }

    // -------------------------------------------------------------------------
    // Field statistics
    // -------------------------------------------------------------------------

    /**
     * Counts grid cells inside the mask that are shallower than fishMaxDepth.
     */
    public static int countAvailableCellsForOxygenUptake(int[] cageDims, double dz, double fishMaxDepth,
                                                         boolean[][][] mask) {
        int res = 0;
        for (int i = 0; i < cageDims[0]; i++)
            for (int j = 0; j < cageDims[1]; j++)
                for (int k = 0; k < cageDims[2]; k++) {
                    double lDepth = ((double) k + 0.5) * dz;
                    if ((mask == null || mask[i][j][k]) && (lDepth < fishMaxDepth)) {
                        res++;
                    }
                }
        return res;
    }

    /**
     * Returns {minValue, meanValue, fractionHypoxic} over masked cells.
     * Hypoxia threshold is {@link #HYPOXIA_THRESHOLD}.
     */
    public static double[] minValueMeanAndFracHypoxia(double[][][] h, boolean[][][] mask) {
        double minval = Double.MAX_VALUE;
        double meanVal = 0;
        int hypoCells = 0, nonHypoCells = 0, cageCells = 0;
        for (int i = 0; i < h.length; i++)
            for (int j = 0; j < h[i].length; j++)
                for (int k = 0; k < h[i][j].length; k++) {
                    if (mask[i][j][k]) {
                        cageCells++;
                        meanVal += h[i][j][k];
                        if (h[i][j][k] < minval)
                            minval = h[i][j][k];
                        if (h[i][j][k] < HYPOXIA_THRESHOLD)
                            hypoCells++;
                        else
                            nonHypoCells++;
                    }
                }
        return new double[]{minval, meanVal / (double) cageCells,
                (double) hypoCells / (double) (nonHypoCells + hypoCells)};
    }

    /**
     * Returns the standard deviation of the 3×3 horizontal neighbourhood around pos,
     * at the same vertical level (pos[2]).
     */
    public static double getStdAround(double[][][] field, int[] pos) {
        double[] values = new double[9];
        values[0] = field[pos[0] - 1][pos[1] - 1][pos[2]];
        values[1] = field[pos[0] - 1][pos[1]][pos[2]];
        values[2] = field[pos[0] - 1][pos[1] + 1][pos[2]];
        values[3] = field[pos[0]][pos[1] - 1][pos[2]];
        values[4] = field[pos[0]][pos[1]][pos[2]];
        values[5] = field[pos[0]][pos[1] + 1][pos[2]];
        values[6] = field[pos[0] + 1][pos[1] - 1][pos[2]];
        values[7] = field[pos[0] + 1][pos[1]][pos[2]];
        values[8] = field[pos[0] + 1][pos[1] + 1][pos[2]];
        return standardDev(values);
    }

    /**
     * Sample standard deviation of a data array.
     */
    public static double standardDev(double[] data) {
        double mean = 0.0;
        for (double v : data) mean += v;
        mean /= data.length;
        double variance = 0.0;
        for (double v : data) variance += Math.pow(v - mean, 2);
        variance /= (data.length - 1);
        return Math.sqrt(variance);
    }

    // -------------------------------------------------------------------------
    // Array utilities
    // -------------------------------------------------------------------------

    /**
     * Fills every cell of a 4-D array with the values from the 1-D array {@code v}.
     * The innermost dimension of {@code h} must be the same length as {@code v}.
     */
    public static void setAllVal(double[][][][] h, double[] v) {
        for (int i = 0; i < h.length; i++)
            for (int j = 0; j < h[i].length; j++)
                for (int k = 0; k < h[i][j].length; k++)
                    for (int z = 0; z < v.length; z++)
                        h[i][j][k][z] = v[z];
    }
}
